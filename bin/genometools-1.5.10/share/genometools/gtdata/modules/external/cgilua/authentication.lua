-- CGILua authentication Module
-- Author: Leonardo Godinho
--
-- Offers a basic API for authentication assuming the presence of 
-- cgilua.POST.user, cgilua.POST.pass, cgilua.QUERY.logout, cgilua.QUERY[tokenName]
--
-- The API consists in the functions
--   check(username, passwd) - Checks if the pair username/passwd is authenticated by the configured method
--   checkURL() - returns the URL for the checking script
--   configure(options, methods) - configures the authentication framework (see /examples/authentication_conf.lua)
--   logoutURL() - returns the URL for the logout page
--   refURL() - returns the original URL being checked for authentication
--   username() - returns the authenticated user if existent
--
-- The authenticated user can be persisted by using a cookie or an ID in the URL
--
-- $Id: authentication.lua,v 1.2 2007/12/05 19:41:13 carregal Exp $


local mime=require"mime" -- from LuaSocket
local md5=require"md5"
require"cgilua.cookies"

local cgilua = cgilua
local string = string
local math = math
local error = error

module"cgilua.authentication"

local authenticatedUser
local configuration
local _check -- check function provided by the configuration

-- Callback functions to manipulate Tokens in the URL
-- if not defined, CGILua standard URLs are assumed

-- Returns the current token
function getToken()
    return cgilua.QUERY[configuration.tokenName]
end

-- Sets the current token
function setToken(token)
    cgilua.QUERY[configuration.tokenName] = token
end
    
-- Returns the current URL
function currentURL()
    local script_name = cgilua.servervariable"SCRIPT_NAME"
    local path_info = cgilua.servervariable"PATH_INFO" or ""
    local query_string = cgilua.servervariable"QUERY_STRING" or ""
    if query_string ~= "" then
        query_string = "?"..query_string
    end
    return cgilua.mkabsoluteurl(script_name..path_info..query_string)
end

-- URL Base64 encoder and decoder functions
-- (http://en.wikipedia.org/wiki/Base64)
--
-- '=' is replaced by ''
-- '+' and '/' are respectively replaced by '*' and '-'
function encodeURLbase64(str)
	local b64str = mime.b64(str)
	local urlb64str = string.gsub(b64str,"=","")
	urlb64str = string.gsub(urlb64str,"+","*")
	urlb64str = string.gsub(urlb64str,"/","-")
	urlb64str = string.gsub(urlb64str," ","_")
	return urlb64str
end

function decodeURLbase64(urlb64str)
	local b64str = string.gsub(urlb64str,"*","+")
	b64str = string.gsub(b64str,"-","/")
	b64str = string.gsub(b64str,"_"," ")
	local b64strPadLen = math.fmod(4 - math.fmod(string.len(b64str), 4), 4)
	b64str = b64str..string.rep("=", b64strPadLen)
	local str = mime.unb64(b64str)
	return str
end

-- Returns the authenticated username or nil if no user is authenticated
function username()
	if authenticatedUser == nil then
		local authenticatedUserData
        local token
		if configuration.tokenPersistence == "url" then
            token = getToken()
		elseif configuration.tokenPersistence == "cookie" then
			token = cgilua.cookies.get(configuration.tokenName)
		end
        if token then
            authenticatedUserData = md5.decrypt(decodeURLbase64(token), configuration.criptKey)
            -- check if IP in crypted data match with client IP
            local authenticatedUserIP = authenticatedUserData and string.gsub(authenticatedUserData, ",.*$","") or nil
            if authenticatedUserIP ~= cgilua.servervariable("REMOTE_ADDR") then
                return nil
            end
            authenticatedUser=authenticatedUserData and string.gsub(authenticatedUserData, "^.*,", "") or nil
        end
	end
	return authenticatedUser
end

-- encrypt the user IP and username for the user hash token
local function cryptUserData()
    if authenticatedUser then
        local userData = cgilua.servervariable("REMOTE_ADDR") ..",".. authenticatedUser
        local cryptedUserData = encodeURLbase64(md5.crypt(userData, configuration.criptKey))
        return cryptedUserData
    end
end

-- defines the logged user name and sets the user hash token
local function setUser(username)
	authenticatedUser = username
    if username then
        local cryptedUserData = cryptUserData()
        if configuration.tokenPersistence == "url" then
            setToken(cryptedUserData)
            cgilua.cookies.delete(configuration.tokenName) -- removes an eventual previous cookie token
        elseif configuration.tokenPersistence == "cookie" then
            cgilua.cookies.set(configuration.tokenName, cryptedUserData)
            setToken() -- remove an eventual previous token from the URLs
        end
    end
end

-- User logout, clear everything
function logout()
    setUser()
    cgilua.cookies.delete(configuration.tokenName)
    setToken()
    cgilua.QUERY.logout = nil
end

-- Checks if a user name/password is authenticated by the configured method
-- if the user is authenticaded then login the user else logout the user
-- returns true if the user has been succesfully authenticated or false plus
-- an error message when the authentication fails
function check(name, pass)
    name = name or cgilua.POST.user
    pass = pass or cgilua.POST.pass
	if name then
        -- Tries to authenticate the user using the configured method
		local retauth,errauth = _check(name, pass)
		if  retauth then
			setUser(name)
			return true
		else
			logout()
			return false, errauth
		end
	else
		local authuser = username()
		if authuser then
			if cgilua.QUERY.logout ~= nil then
				logout()
				return false
			end
		end
		return authuser
	end
end

-- Returns a authentication URL with ref URL as a parameter,
-- accepts an optional value for the logout action
function checkURL(ref, tologout)
    local token
    if configuration.tokenPersistence == "url" then
        token = getToken()
    elseif configuration.tokenPersistence == "cookie" then
        token = cgilua.cookies.get(configuration.tokenName)
    end

    -- As HTTP header referer information can violate privacy, 
    -- some browsers allow the user to disable the sending of referer information.
    -- Some proxy and firewall software will also filter out referer information,
    -- to avoid leaking the location of non-public websites.
    -- So we send the current URL as an URL parameter to the login URL.     
	setToken()
	local args = {ref = ref or currentURL(), logout = tologout}
	if string.find(configuration.checkURL, "^https?:") then
		local params = "?"..urlcode.encodetable(args)
		return configuration.checkURL..params
	end
	return cgilua.mkabsoluteurl(cgilua.mkurlpath(configuration.checkURL, args))
end

-- Returns the logout URL, based on the login URL
function logoutURL()
	return checkURL(nil, 1)
end

-- Returns the referenced URL, the one supposed to be offered only for authenticated users
function refURL()
    local url
    local baseURL = cgilua.QUERY.ref or configuration.checkURL
    if string.find(baseURL, "\?") then
        url = string.gsub(baseURL, "\?", "?"..configuration.tokenName.."="..cryptUserData().."&")
    else
        url = baseURL.."?"..configuration.tokenName.."="..cryptUserData()
    end
	return url
end

-- Sets the current configuration
function configure(options, methods)
    configuration = options
    local method = methods[options.method] or {}
    
    if method.check then
        _check = method.check
    end
    if method.username then
        username = method.username
    end
end