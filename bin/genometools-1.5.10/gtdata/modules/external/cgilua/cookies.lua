----------------------------------------------------------------------------
-- Cookies Library
--
-- @release $Id: cookies.lua,v 1.8 2008/04/24 13:42:04 mascarenhas Exp $
----------------------------------------------------------------------------

require"cgilua.urlcode"

local error = error
local format, gsub, strfind = string.format, string.gsub, string.find
local date = os.date
local escape, unescape = cgilua.urlcode.escape, cgilua.urlcode.unescape
local function header(...)
   return SAPI.Response.header(...)
end
local function write(...)
   return SAPI.Response.write(...)
end
local function servervariable(...)
   return SAPI.Request.servervariable(...)
end

module ("cgilua.cookies")

local function optional (what, name)
  if name ~= nil and name ~= "" then
    return format("; %s=%s", what, name)
  else
    return ""
  end
end


local function build (name, value, options)
  if not name or not value then
    error("cookie needs a name and a value")
  end
  local cookie = name .. "=" .. escape(value)
  options = options or {}
  if options.expires then
    local t = date("!%A, %d-%b-%Y %H:%M:%S GMT", options.expires)
    cookie = cookie .. optional("expires", t)
  end
  cookie = cookie .. optional("path", options.path)
  cookie = cookie .. optional("domain", options.domain)
  cookie = cookie .. optional("secure", options.secure)
  return cookie
end


----------------------------------------------------------------------------
-- Sets a value to a cookie, with the given options.
-- Generates a header "Set-Cookie", thus it can only be used in Lua Scripts.
-- @param name String with the name of the cookie.
-- @param value String with the value of the cookie.
-- @param options Table with the options (optional).

function set (name, value, options)
  header("Set-Cookie", build(name, value, options))
end


----------------------------------------------------------------------------
-- Sets a value to a cookie, with the given options.
-- Generates an HTML META tag, thus it can be used in Lua Pages.
-- @param name String with the name of the cookie.
-- @param value String with the value of the cookie.
-- @param options Table with the options (optional).

function sethtml (name, value, options)
  write(format('<meta http-equiv="Set-Cookie" content="%s">', 
                build(name, value, options)))
end


----------------------------------------------------------------------------
-- Gets the value of a cookie.
-- @param name String with the name of the cookie.
-- @return String with the value associated with the cookie.

function get (name)
  local cookies = servervariable"HTTP_COOKIE" or ""
  cookies = ";" .. cookies .. ";"
  cookies = gsub(cookies, "%s*;%s*", ";")   -- remove extra spaces
  local pattern = ";" .. name .. "=(.-);"
  local _, __, value = strfind(cookies, pattern)
  return value and unescape(value)
end


----------------------------------------------------------------------------
-- Deletes a cookie, by setting its value to "xxx".
-- @param name String with the name of the cookie.
-- @param options Table with the options (optional).

function delete (name, options)
  options = options or {}
  options.expires = 1
  set(name, "xxx", options)
end
