-- CGILua dispatcher module
-- @release $Id: dispatcher.lua,v 1.8 2007/12/07 18:49:49 carregal Exp $

module(..., package.seeall)

-- Checks if an URL matches a route pattern
local function route_match(url, pattern) 
    local params = {}
    local captures = string.gsub(pattern, "(/$[%w_-]+)", "/([^/]*)")
    local url_parts = {string.match(url, captures)}
    local i = 1
    for name in string.gmatch(pattern, "/$([%w_-]+)") do
        params[name] = url_parts[i]
        i = i + 1
    end
    return next(params) and params
end

local route_URLs = {}

-- Maps the correct function for an URL
local function route_map(url) 
    for i, v in ipairs(route_URLs) do
        local pattern, f, name = unpack(v)
        local params = route_match(url, pattern)
        if params then 
            return f, params 
        end
    end
end

-- Returns an URL for a named route
-- @param map_name Name associated with the map in the routed URL table.
-- @param params Table of named parameters used in the URL map
-- @param queryargs Optional table of named parameters used for the QUERY part of the URL
function route_url(map_name, params, queryargs)
	local queryparams = ""
	if queryargs then
		queryparams = "?"..cgilua.urlcode.encodetable(queryargs)
	end
	for i, v in ipairs(route_URLs) do
        local pattern, f, name = unpack(v)
        if name == map_name then
            local url = string.gsub(pattern, "$([%w_-]+)", params)
            url = cgilua.urlpath.."/"..cgilua.app_name..url..queryparams
            return url
        end
    end
end

-- Defines the routing using a table of URLs maps or a single map
-- a map defines a URL mask using $name to extract parameters,
-- a function to be called with the extracted parameters and
-- a name for the map when used with route_url
-- @param table of maps or a single map
function route(URLs)
	URLs = URLs or {}
	if type(URLs[1]) == "string" then
		-- accepts a single map as the only entry in a map table
		URLs = {URLs}
	end
    route_URLs = URLs
    f, args = route_map(cgilua.script_vpath)

    if f then
        return f(args)
    else
        error("Missing page parameters")
    end
end