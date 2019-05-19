----------------------------------------------------------------------------
-- Process POST data.
-- This library depends on some functions that read POST data and other
-- HTTP information.  A beginning is:
--	require"post"
--	local params = {}
--	post.parsedata {
--		read = ap.get_client_block or io.input,
--		discardinput = ap.discard_request_body,
--		content_type = ap.get_header"content-type" or os.getenv"CONTENT_TYPE",
--		content_length = ap.get_header"content-length" or os.getenv"CONTENT_LENGTH",
--		maxinput = 1024 * 1024,
--		maxfilesize = 512 * 1024,
--		args = params,
--	}
--
-- @release $Id: post.lua,v 1.17 2008/04/03 21:55:28 mascarenhas Exp $
----------------------------------------------------------------------------

require"cgilua.readuntil"
require"cgilua.urlcode"

local assert, error, pairs, tonumber, tostring, type = assert, error, pairs, tonumber, tostring, type
local getn, tinsert = table.getn, table.insert
local format, gsub, strfind, strlower, strlen = string.format, string.gsub, string.find, string.lower, string.len
local min = math.min
local iterate = cgilua.readuntil.iterate
local urlcode = cgilua.urlcode
local tmpfile = cgilua.tmpfile

-- environment for processing multipart/form-data input
local boundary = nil      -- boundary string that separates each 'part' of input
local maxfilesize = nil   -- maximum size for file upload
local maxinput = nil      -- maximum size of total POST data
local inputfile = nil     -- temporary file for inputting form-data
local bytesleft = nil     -- number of bytes yet to be read
local content_type = nil  -- request's content-type
-- local functions
local discardinput = nil  -- discard all remaining input
local readuntil = nil     -- read until delimiter
local read = nil          -- basic read function

module ("cgilua.post")

--
-- Extract the boundary string from CONTENT_TYPE metavariable
--
local function getboundary ()
  local _,_,boundary = strfind (content_type, "boundary%=(.-)$")
  return  "--"..boundary 
end

--
-- Create a table containing the headers of a multipart/form-data field
--
local function breakheaders (hdrdata)
  local headers = {}
  gsub (hdrdata, '([^%c%s:]+):%s+([^\n]+)', function(type,val)
    type = strlower(type)
    headers[type] = val
  end)
  return headers
end

--
-- Read the headers of the next multipart/form-data field 
--
--  This function returns a table containing the headers values. Each header
--  value is indexed by the corresponding header "type". 
--  If end of input is reached (no more fields to process) it returns nil.
--
local function readfieldheaders ()
	local EOH = "\r\n\r\n" -- <CR><LF><CR><LF>
	local hdrdata = ""
	local out = function (str) hdrdata = hdrdata..str end
	if readuntil (EOH, out) then
		-- parse headers
		return breakheaders (hdrdata)
	else
		-- no header found
		return nil
	end
end

--
-- Extract a field name (and possible filename) from its disposition header
--
local function getfieldnames (headers)
  local disposition_hdr = headers["content-disposition"]
  local attrs = {}
  if disposition_hdr then
    gsub(disposition_hdr, ';%s*([^%s=]+)="(.-)"', function(attr, val)
	   attrs[attr] = val
         end)
  else
    error("Error processing multipart/form-data."..
          "\nMissing content-disposition header")
  end
  return attrs.name, attrs.filename
end

--
-- Read the contents of a 'regular' field to a string
--
local function readfieldcontents ()
	local value = ""
	local boundaryline = "\r\n"..boundary
	local out = function (str) value = value..str end
	if readuntil (boundaryline, out) then
		return value
	else
		error("Error processing multipart/form-data.\nUnexpected end of input\n")
	end
end

--
-- Read the contents of a 'file' field to a temporary file (file upload)
--
local function fileupload (filename)
	-- create a temporary file for uploading the file field
	local file, err = tmpfile()
	if file == nil then
		discardinput(bytesleft)
		error("Cannot create a temporary file.\n"..err)
	end      
	local bytesread = 0
	local boundaryline = "\r\n"..boundary
	local out = function (str)
		local sl = strlen (str)
		if bytesread + sl > maxfilesize then
			discardinput (bytesleft)
			error (format ("Maximum file size (%d kbytes) exceeded while uploading `%s'", maxfilesize / 1024, filename))
		end
		file:write (str)
		bytesread = bytesread + sl
	end
	if readuntil (boundaryline, out) then
		file:seek ("set", 0)
		return file, bytesread
	else
		error (format ("Error processing multipart/form-data.\nUnexpected end of input while uploading %s", filename))
	end
end

--
-- Compose a file field 'value' 
--
local function filevalue (filehandle, filename, filesize, headers)
  -- the temporary file handle
  local value = { file = filehandle,
                  filename = filename,
                  filesize = filesize }
  -- copy additional header values
  for hdr, hdrval in pairs(headers) do
    if hdr ~= "content-disposition" then
      value[hdr] = hdrval
    end
  end
  return value
end

--
-- Process multipart/form-data 
--
-- This function receives the total size of the incoming multipart/form-data, 
-- the maximum size for a file upload, and a reference to a table where the 
-- form fields should be stored.
--
-- For every field in the incoming form-data a (name=value) pair is 
-- inserted into the given table. [[name]] is the field name extracted
-- from the content-disposition header.
--
-- If a field is of type 'file' (i.e., a 'filename' attribute was found
-- in its content-disposition header) a temporary file is created 
-- and the field contents are written to it. In this case,
-- [[value]] has a table that contains the temporary file handle 
-- (key 'file') and the file name (key 'filename'). Optional headers
-- included in the field description are also inserted into this table,
-- as (header_type=value) pairs.
--
-- If the field is not of type 'file', [[value]] contains the field 
-- contents.
--
local function Main (inputsize, args)

	-- create a temporary file for processing input data
	local inputf,err = tmpfile()
	if inputf == nil then
		discardinput(inputsize)
		error("Cannot create a temporary file.\n"..err)
	end

	-- set the environment for processing the multipart/form-data
	inputfile = inputf
	bytesleft = inputsize
	maxfilesize = maxfilesize or inputsize 
	boundary = getboundary()

	while true do
		-- read the next field header(s)
		local headers = readfieldheaders()
		if not headers then break end	-- end of input

		-- get the name attributes for the form field (name and filename)
		local name, filename = getfieldnames(headers)

		-- get the field contents
		local value
		if filename then
			local filehandle, filesize = fileupload(filename)
			value = filevalue(filehandle, filename, filesize, headers)
		else
			value = readfieldcontents()
		end

		-- insert the form field into table [[args]]
		urlcode.insertfield(args, name, value)
	end
end

--
-- Initialize the library by setting the dependent functions:
--	content_type            = value of "Content-type" header
--	content_length          = value of "Content-length" header
--	read                    = function that can read POST data
--	discardinput (optional) = function that discard POST data
--	maxinput (optional)     = limit of POST data (in bytes)
--	maxfilesize (optional)  = limit of uploaded file(s) (in bytes)
--
local function init (defs)
	assert (defs.read)
	read = defs.read
	readuntil = iterate (function ()
        if bytesleft then
            if bytesleft <= 0 then return nil end
            local n = min (bytesleft, 2^13) -- 2^13 == 8192
            local bytes = read (n)
            bytesleft = bytesleft - #bytes
            return bytes
        end
	end)
	if defs.discard_function then
		discardinput = defs.discardinput
	else
		discardinput = function (inputsize)
			readuntil ('\0', function()end)
		end
	end
	content_type = defs.content_type
	if defs.maxinput then
		maxinput = defs.maxinput
	end
	if defs.maxfilesize then
		maxfilesize = defs.maxfilesize
	end
end

----------------------------------------------------------------------------
-- Parse the POST REQUEST incoming data according to its "content type"
-- as defined by the metavariable CONTENT_TYPE (RFC CGI)
--
--  An error is issued if the "total" size of the incoming data
--   (defined by the metavariable CONTENT_LENGTH) exceeds the
--   maximum input size allowed
----------------------------------------------------------------------------
function parsedata (defs)
	assert (type(defs.args) == "table", "field `args' must be a table")
	init (defs)
	-- get the "total" size of the incoming data
	local inputsize = tonumber(defs.content_length) or 0
	if inputsize > maxinput then
		-- some Web Servers (like IIS) require that all the incoming data is read 
		bytesleft = inputsize
		discardinput(inputsize)
		error(format("Total size of incoming data (%d KB) exceeds configured maximum (%d KB)",
			inputsize /1024, maxinput / 1024))
	end

	-- process the incoming data according to its content type
	local contenttype = content_type
	if not contenttype then
		error("Undefined Media Type") 
	end
	if strfind(contenttype, "x-www-form-urlencoded", 1, true) then
		urlcode.parsequery (read (inputsize), defs.args)
	elseif strfind(contenttype, "multipart/form-data", 1, true) then
		Main (inputsize, defs.args)
	elseif strfind (contenttype, "application/xml", 1, true) or strfind (contenttype, "text/xml", 1, true) or strfind (contenttype, "text/plain", 1, true) then
		tinsert (defs.args, read (inputsize))
	else
		error("Unsupported Media Type: "..contenttype)
	end
end
