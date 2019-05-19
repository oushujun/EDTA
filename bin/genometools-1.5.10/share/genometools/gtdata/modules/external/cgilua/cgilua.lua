----------------------------------------------------------------------------
-- CGILua library.
--
-- @release $Id: cgilua.lua,v 1.81 2008/05/19 18:13:36 carregal Exp $
----------------------------------------------------------------------------

local _G, SAPI = _G, SAPI
local urlcode = require"cgilua.urlcode"
local lp = require"cgilua.lp"
local lfs = require"lfs"
local debug = require"debug"
local assert, error, ipairs, select, tostring, type, unpack, xpcall = assert, error, ipairs, select, tostring, type, unpack, xpcall
local pairs = pairs
local gsub, format, strfind, strlower, strsub, match = string.gsub, string.format, string.find, string.lower, string.sub, string.match
local setmetatable = setmetatable
local _open = io.open
local tinsert, tremove, concat = table.insert, table.remove, table.concat
local foreachi = table.foreachi
local date = os.date
local os_tmpname = os.tmpname
local getenv = os.getenv
local remove = os.remove
local seeall = package.seeall
local setfenv = setfenv

lp.setoutfunc ("cgilua.put")
lp.setcompatmode (true)

module ("cgilua")

_COPYRIGHT = "Copyright (C) 2003 Kepler Project"
_DESCRIPTION = "CGILua is a tool for creating dynamic Web pages and manipulating input data from forms"
_VERSION = "CGILua 5.1.2"
--
-- Internal state variables.
local _default_errorhandler = debug.traceback
local _errorhandler = _default_errorhandler
local _default_erroroutput = function (msg)

    if type(msg) ~= "string" and type(msg) ~= "number" then
        msg = format ("bad argument #1 to 'error' (string expected, got %s)", type(msg))
    end
  
	-- Logging error
	SAPI.Response.errorlog (msg)
	SAPI.Response.errorlog (" ")

	SAPI.Response.errorlog (SAPI.Request.servervariable"REMOTE_ADDR")
	SAPI.Response.errorlog (" ")

	SAPI.Response.errorlog (date())
	SAPI.Response.errorlog ("\n")

	-- Building user message
	msg = gsub (gsub (msg, "\n", "<br>\n"), "\t", "&nbsp;&nbsp;")
	SAPI.Response.contenttype ("text/html")
	SAPI.Response.write ("<html><head><title>CGILua Error</title></head><body>" .. msg .. "</body></html>")
end
local _erroroutput = _default_erroroutput
local _default_maxfilesize = 512 * 1024
local _maxfilesize = _default_maxfilesize
local _default_maxinput = 1024 * 1024
local _maxinput = _default_maxinput
script_path = false

--
-- Header functions

----------------------------------------------------------------------------
-- Sends a header.
-- @param header String with the header.
-- @param value String with the corresponding value.
----------------------------------------------------------------------------
function header(...)
   return SAPI.Response.header(...)
end

----------------------------------------------------------------------------
-- Sends a Content-type header.
-- @param type String with the type of the header.
-- @param subtype String with the subtype of the header.
----------------------------------------------------------------------------
function contentheader (type, subtype)
	SAPI.Response.contenttype (type..'/'..subtype)
end

----------------------------------------------------------------------------
-- Sends the HTTP header "text/html".
----------------------------------------------------------------------------
function htmlheader()
	SAPI.Response.contenttype ("text/html")
end
local htmlheader = htmlheader

----------------------------------------------------------------------------
-- Sends an HTTP header redirecting the browser to another URL
-- @param url String with the URL.
-- @param args Table with the arguments (optional).
----------------------------------------------------------------------------
function redirect (url, args)
	if strfind(url,"^https?:") then
		local params=""
		if args then
			params = "?"..urlcode.encodetable(args)
		end
		return SAPI.Response.redirect(url..params)
	else
		return SAPI.Response.redirect(mkabsoluteurl(mkurlpath(url,args)))
	end
end

----------------------------------------------------------------------------
-- Returns a server variable
-- @param name String with the name of the server variable.
-- @return String with the value of the server variable.
----------------------------------------------------------------------------
function servervariable(...)
   return SAPI.Request.servervariable(...)
end

----------------------------------------------------------------------------
-- Primitive error output function
-- @param msg String (or number) with the message.
-- @param level String with the error level (optional).
----------------------------------------------------------------------------
function errorlog (msg, level)
	local t = type(msg)
	if t == "string" or t == "number" then
		SAPI.Response.errorlog (msg, level)
	else
		error ("bad argument #1 to `cgilua.errorlog' (string expected, got "..t..")", 2)
	end
end

----------------------------------------------------------------------------
-- Converts all its arguments to strings before sending them to the server.
----------------------------------------------------------------------------
function print (...)
	local args = { ... }
	for i = 1, select("#",...) do
		args[i] = tostring(args[i])
	end
	SAPI.Response.write (concat(args,"\t"))
	SAPI.Response.write ("\n")
end

----------------------------------------------------------------------------
-- Function 'put' sends its arguments (basically strings of HTML text)
--  to the server
-- Its basic implementation is to use Lua function 'write', which writes
--  each of its arguments (strings or numbers) to file _OUTPUT (a file
--  handle initialized with the file descriptor for stdout)
-- @param s String (or number) with output.
----------------------------------------------------------------------------
function put (...)
   return SAPI.Response.write(...)
end

-- Returns the current errorhandler
function _geterrorhandler(msg)
    return _errorhandler(msg)
end

--
-- Executes a function using the CGILua error handler.
-- @param f Function to be called.
--
function pcall (f)
	local results = {xpcall (f, _geterrorhandler)}
	local ok = results[1]
	tremove(results, 1)
	if ok then
	        if #results == 0 then results = { true } end
		return unpack(results)
	else
		_erroroutput (unpack(results))
	end
end

local function buildscriptenv()
  local env = { print = _M.print, write = _M.put }
  setmetatable(env, { __index = _G, __newindex = _G })
  return env
end

----------------------------------------------------------------------------
-- Execute a script
--  If an error is found, Lua's error handler is called and this function
--  does not return
-- @param filename String with the name of the file to be processed.
-- @return The result of the execution of the file.
----------------------------------------------------------------------------
function doscript (filename)
  local f, err = _G.loadfile(filename)
  if not f then
    error (format ("Cannot execute `%s'. Exiting.\n%s", filename, err))
  else
    local env = buildscriptenv()
    setfenv(f, env)
    return pcall(f)
  end
end

----------------------------------------------------------------------------
-- Execute the file if there is no "file error".
--  If an error is found, and it is not a "file error", Lua 'error'
--  is called and this function does not return
-- @param filename String with the name of the file to be processed.
-- @return The result of the execution of the file or nil (in case the
--	file does not exists or if it cannot be opened).
-- @return It could return an error message if the file cannot be opened.
----------------------------------------------------------------------------
function doif (filename)
	if not filename then return end    -- no file
	local f, err = _open(filename)
	if not f then return nil, err end    -- no file (or unreadable file)
	f:close()
	return doscript (filename)
end

---------------------------------------------------------------------------
-- Set the maximum "total" input size allowed (in bytes)
-- @param nbytes Number of the maximum size (in bytes) of the whole POST data.
---------------------------------------------------------------------------
function setmaxinput(nbytes)
	_maxinput = nbytes
end

---------------------------------------------------------------------------
-- Set the maximum size for an "uploaded" file (in bytes)
-- Might be less or equal than _maxinput.
-- @param nbytes Number of the maximum size (in bytes) of a file.
---------------------------------------------------------------------------
function setmaxfilesize(nbytes)
	_maxfilesize = nbytes
end


-- Default path for temporary files
tmp_path = _G.CGILUA_TMP or getenv("TEMP") or getenv ("TMP") or "/tmp"

-- Default function for temporary names
-- @returns a temporay name using os.tmpname
tmpname = function()
    local tempname = os_tmpname()
    -- Lua os.tmpname returns a full path in Unix, but not in Windows
    -- so we strip the eventual prefix
    tempname = gsub(tempname, "(/tmp/)", "")
    return tempname
end

local _tmpfiles = {}

---------------------------------------------------------------------------
-- Returns a temporary file in a directory using a name generator
-- @param dir Base directory for the temporary file
-- @param namefunction Name generator function
---------------------------------------------------------------------------
function tmpfile(dir, namefunction)
    dir = dir or tmp_path
    namefunction = namefunction or tmpname
    local tempname = namefunction()
    local filename = dir.."/"..tempname
    local file, err = _open(filename, "wb+")
    if file then
        tinsert(_tmpfiles, {name = filename, file = file})
    end
    return file, err
end


----------------------------------------------------------------------------
-- Preprocess the content of a mixed HTML file and output a complete
--   HTML document ( a 'Content-type' header is inserted before the
--   preprocessed HTML )
-- @param filename String with the name of the file to be processed.
-- @param env Optional environment
----------------------------------------------------------------------------
function handlelp (filename, env)
  env = env or buildscriptenv()
  htmlheader ()
  lp.include (filename, env)
end

----------------------------------------------------------------------------
-- Builds a handler that sends a header and the contents of the given file.
-- Sends the contents of the file to the output without processing it.
-- @param type String with the type of the header.
-- @param subtype String with the subtype of the header.
-- @return Function (which receives a filename as argument) that produces
--	the header and copies the content of the given file.
----------------------------------------------------------------------------
function buildplainhandler (type, subtype)
	return function (filename)
		local fh, err = _open (filename, "rb")
        local contents = ""
        if fh then
            contents = fh:read("*a")
            fh:close()
        else
            error(err)
        end
        header("Content-Lenght", #contents)
		contentheader (type, subtype)
		put (contents)
	end
end

----------------------------------------------------------------------------
-- Builds a handler that sends a header and the processed file.
-- Processes the file as a Lua Page.
-- @param type String with the type of the header.
-- @param subtype String with the subtype of the header.
-- @return Function (which receives a filename as argument) that produces
--	the header and processes the given file.
----------------------------------------------------------------------------
function buildprocesshandler (type, subtype)
	return function (filename)
		 local env = buildscriptenv()
		 contentheader (type, subtype)
		 lp.include (filename, env)
	end
end

----------------------------------------------------------------------------
-- Builds the default handler table from cgilua.mime
----------------------------------------------------------------------------
local function buildhandlers()
    local mime = _G.require "cgilua.mime"
    for ext, mediatype in pairs(mime) do
        local t, st = match(mediatype, "([^/]*)/([^/]*)")
        addscripthandler(ext, buildplainhandler(t, st))
    end
end

----------------------------------------------------------------------------
-- Create an URL path to be used as a link to a CGILua script
-- @param script String with the name of the script.
-- @param args Table with arguments to script (optional).
-- @return String in URL format.
----------------------------------------------------------------------------
function mkurlpath (script, args)
	-- URL-encode the parameters to be passed do the script
	local params = ""
	if args then
		params = "?"..urlcode.encodetable(args)
	end
	if strsub(script,1,1) == "/" then
		return urlpath .. script .. params
	else
		return urlpath .. script_vdir .. script .. params
	end
end

----------------------------------------------------------------------------
-- Create an absolute URL containing the given URL path
-- @param path String with the path.
-- @param protocol String with the name of the protocol (default = "http").
-- @return String in URL format.
----------------------------------------------------------------------------
function mkabsoluteurl (path, protocol)
	protocol = protocol or "http"
	if path:sub(1,1) ~= '/' then
		path = '/'..path
	end
	return format("%s://%s:%s%s",
		protocol,
		servervariable"SERVER_NAME",
		servervariable"SERVER_PORT",
		path)
end

----------------------------------------------------------------------------
-- Extract the "directory" and "file" parts of a path
-- @param path String with a path.
-- @return String with the directory part.
-- @return String with the file part.
----------------------------------------------------------------------------
function splitonlast (path, sep)
	local dir,file = match(path,"^(.-)([^:/\\]*)$")
	return dir,file
end

splitpath = splitonlast -- compatibility with previous versions

----------------------------------------------------------------------------
-- Extracts the first and remaining parts of a path
-- @param path separator (defaults to "/")
-- @return String with the extracted part.
-- @return String with the remaining path.
----------------------------------------------------------------------------
function splitonfirst(path, sep)
    local first, rest = match(path, "^/([^:/\\]*)(.*)")
    return first, rest
end

--
-- Define variables and build the cgilua.POST, cgilua.GET tables.
--
local function getparams ()
    requestmethod = servervariable"REQUEST_METHOD"
	-- Fill in the POST table.
	POST = {}
	if  requestmethod == "POST" then
		post.parsedata {
			read = SAPI.Request.getpostdata,
			discardinput = ap and ap.discard_request_body,
			content_type = servervariable"CONTENT_TYPE",
			content_length = servervariable"CONTENT_LENGTH",
			maxinput = _maxinput,
			maxfilesize = _maxfilesize,
			args = POST,
		}
	end
	-- Fill in the QUERY table.
	QUERY = {}
	urlcode.parsequery (servervariable"QUERY_STRING", QUERY)
end

--
-- Stores all script handlers and the file extensions used to identify
-- them. Loads the default 
local _script_handlers = {}
--
-- Default handler.
-- Sends the contents of the file to the output without processing it.
-- This relies in the browser being able to discover the content type
-- which is not reliable.
-- @param filename String with the name of the file.
--
local function default_handler (filename)
	local fh, err = _open (filename, "rb")
    local contents
    if fh then
        contents = fh:read("*a")
        fh:close()
    else
        error(err)
    end
    header("Content-Lenght", #contents)
    put ("\n")
	put (contents)
end

----------------------------------------------------------------------------
-- Add a script handler.
-- @param file_extension String with the lower-case extension of the script.
-- @param func Function to handle this kind of scripts.
----------------------------------------------------------------------------
function addscripthandler (file_extension, func)
	assert (type(file_extension) == "string", "File extension must be a string")
	if strfind (file_extension, '%.', 1) then
		file_extension = strsub (file_extension, 2)
	end
	file_extension = strlower(file_extension)
	assert (type(func) == "function", "Handler must be a function")

	_script_handlers[file_extension] = func
end

---------------------------------------------------------------------------
-- Obtains the handler corresponding to the given script path.
-- @param path String with a script path.
-- @return Function that handles it or nil.
----------------------------------------------------------------------------
function getscripthandler (path)
	local i,f, ext = strfind (path, "%.([^.]+)$")
	return _script_handlers[strlower(ext or '')]
end

---------------------------------------------------------------------------
-- Execute the given path with the corresponding handler.
-- @param path String with a script path.
-- @return The returned values from the script.
---------------------------------------------------------------------------
function handle (path)
	local h = getscripthandler (path) or default_handler
    return h (path)
end

---------------------------------------------------------------------------
-- Sets "errorhandler" function
-- This function is called by Lua when an error occurs.
-- It receives the error message generated by Lua and it is resposible
-- for the final message which should be returned.
-- @param Function.
---------------------------------------------------------------------------
function seterrorhandler (f)
	local tf = type(f)
	if tf == "function" then
		_errorhandler = f
	else
		error (format ("Invalid type: expected `function', got `%s'", tf))
	end
end

---------------------------------------------------------------------------
-- Defines the "erroroutput" function
-- This function is called to generate the error output.
-- @param Function.
---------------------------------------------------------------------------
function seterroroutput (f)
	local tf = type(f)
	if tf == "function" then
		_erroroutput = f
	else
		error (format ("Invalid type: expected `function', got `%s'", tf))
	end
end

--
-- Stores all close functions in order they are set.
local _close_functions = {
}

---------------------------------------------------------------------------
-- Adds a function to be executed after the script.
-- @param f Function to be registered.
---------------------------------------------------------------------------
function addclosefunction (f)
	local tf = type(f)
	if tf == "function" then
		tinsert (_close_functions, f)
	else
		error (format ("Invalid type: expected `function', got `%s'", tf))
	end
end

--
-- Close function.
--
local function close()
	for i = #_close_functions, 1, -1 do
		_close_functions[i]()
	end
end

--
-- Stores all open functions in order they are set.
local _open_functions = {
}

---------------------------------------------------------------------------
-- Adds a function to be executed before the script.
-- @param f Function to be registered.
---------------------------------------------------------------------------
function addopenfunction (f)
	local tf = type(f)
	if tf == "function" then
		tinsert (_open_functions, f)
	else
		error (format ("Invalid type: expected `function', got `%s'", tf))
	end
end

--
-- Open function.
-- Call all defined open-functions in the order they were created.
--
local function open()
	for i = #_open_functions, 1, -1 do
		_open_functions[i]()
	end
end

--
-- Resets CGILua's state.
--
local function reset ()
	script_path = false
	_maxfilesize = _default_maxfilesize
	_maxinput = _default_maxinput
	-- Error Handling
	_errorhandler = _default_errorhandler
	_erroroutput = _default_erroroutput
	-- Handlers
	_script_handlers = {}
	_open_functions = {}
	_close_functions = {}
    -- clean temporary files
    foreachi(_tmpfiles, function (i, v)
        v.file:close()
        local _, err = remove(v.name)
        if err then
            error(err)
        end
    end)
end

---------------------------------------------------------------------------
-- Request processing.
---------------------------------------------------------------------------
function main ()
        SAPI = _G.SAPI
	buildhandlers()    
	-- Default handler values
	addscripthandler ("lua", doscript)
	addscripthandler ("lp", handlelp)
	-- Looks for an optional loader module
	pcall (function () _G.require"cgilua.loader" end)

	-- post.lua needs to be loaded after cgilua.lua is compiled
	pcall (function () _G.require"cgilua.post" end)

	if loader then
		loader.init()
	end
    
	-- Build QUERY/POST tables
	if not pcall (getparams) then return nil end

	local result
	-- Executes the optional loader module
	if loader then
		loader.run()
	end

	-- Changing curent directory to the script's "physical" dir
	local curr_dir = lfs.currentdir ()
	pcall (function () lfs.chdir (script_pdir) end)

	-- Opening functions
	pcall (open)

	-- Executes the script
	result, err = pcall (function () return handle (script_file) end)
    
	-- Closing functions
	pcall (close)
	-- Changing to original directory
	pcall (function () lfs.chdir (curr_dir) end)

	-- Cleanup
	reset ()
	if result then -- script executed ok!
		return result
	end
end
