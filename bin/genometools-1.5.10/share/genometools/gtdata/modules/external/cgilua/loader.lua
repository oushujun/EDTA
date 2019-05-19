-- CGILua loader
--
-- Assumes two global variables pointing to directories:
--  CGILUA_CONF - CGILua configuration directory
--  CGILUA_TMP  - CGILua temporary files directory
--
-- # $Id: loader.lua,v 1.7 2008/03/28 22:01:18 ignacio Exp $

require "lfs"

module("cgilua.loader", package.seeall)

function init()
   -- Loads the user configuration
   cgilua.doif (CGILUA_CONF.."/config.lua")
end

function run()
   local servervariable = cgilua.servervariable
   cgilua.script_vpath = cgilua.script_vpath or servervariable"PATH_INFO"

   if cgilua.script_vpath == nil or cgilua.script_vpath == "" then
      cgilua.script_vpath = "/"
   end

   local document_root = cgilua.pdir or servervariable"DOCUMENT_ROOT"

   if not cgilua.script_path then
      if CGILUA_ISDIRECT then -- the launcher is associated with scripts (Xavante, mod_lua or isapi usually)
	 if servervariable"PATH_TRANSLATED" ~= "" then
	    cgilua.script_path = servervariable"PATH_TRANSLATED"
	 else
	    cgilua.script_path = servervariable"SCRIPT_FILENAME"
	 end
      else -- the launcher executable is present in the URL (CGI or FastCGI only)
	 if document_root == nil or document_root == "" then
	    -- tries to obtain documents_root from the server
	    local path_info = cgilua.script_vpath
	    if path_info == nil or path_info == "" or path_info == "/" then
	       document_root = cgilua.pdir or servervariable"PATH_TRANSLATED"
	    else
	       if string.find(servervariable("SERVER_SOFTWARE"), "IIS") then
		  path_info = string.gsub(path_info, "/", "\\")
	       end
	       document_root = cgilua.pdir or string.gsub(servervariable"PATH_TRANSLATED", path_info, "")
	    end
	 end
	 if cgilua.use_executable_name then
	    -- looks for a Lua script with the same name as the executable
	    local _, name = cgilua.splitpath(servervariable"SCRIPT_NAME")
	    name = string.gsub(name, "%.[^%.]-$","")
		if name and lfs.attributes(document_root.."/"..name..".lua") then
			cgilua.script_path = document_root.."/"..name..".lua"
		end
	 else
	    -- uses /index.lua then /index.lp as the default script
	    if cgilua.script_vpath == "/" then
	       if lfs.attributes(document_root.."/index.lua") then
		  cgilua.script_vpath = "/index.lua"
	       elseif lfs.attributes(document_root.."/index.lp") then
		  cgilua.script_vpath = "/index.lp"
	       else
		  error("Kepler is correctly configured, but you didn't provide a script!")
	       end                    
	    end
	    -- checks if PATH_INFO refers to a valid file and ajusts the settings accordingly
	    local filepath, path_info = string.match (cgilua.script_vpath, "^([^%.]-%.[^/]+)(.*)")
	    if filepath and lfs.attributes(document_root..filepath) then
	       -- if one is found use it
	       cgilua.script_path = document_root..filepath
	       cgilua.script_vpath = path_info
	       cgilua.urlpath = cgilua.urlpath or servervariable"SCRIPT_NAME"..filepath
	    else
	       -- otherwise go with the current PATH_INFO
	       cgilua.script_path = document_root..cgilua.script_vpath
	    end
	 end
      end
   end
   
   -- define other cgilua vars so mkurlpath can work correctly
   if cgilua.script_vpath then
      cgilua.script_vdir = cgilua.splitpath (cgilua.script_vpath)
      cgilua.urlpath = cgilua.urlpath or servervariable"SCRIPT_NAME"
   else
      cgilua.script_vdir = cgilua.splitpath (servervariable"SCRIPT_NAME")
      cgilua.urlpath = cgilua.urlpath or ""
   end

   cgilua.script_pdir, cgilua.script_file = cgilua.splitpath (cgilua.script_path)
end
