function table.contains(tab, element)
  for _, value in pairs(tab) do
    if value == element then
      return true
    end
  end
  return false
end

function table.pretty_array(tab)
  return  "[" .. (table.concat(tab, ", "))  .. "]"
end

function string.char_count(str, char)
    if not str then return 0 end
    local count = 0
    local byte_char = string.byte(char)
    for i = 1, #str do
      if string.byte(str, i) == byte_char then
        count = count + 1
      end
    end
    return count
end

function string.split(str, pat)
   local t = {}
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
        table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end

function collect(iterator)
  local t = {}
  for v in iterator do
    t[#t+1] = v
  end
  return t
end

function count(iterator)
  local t = 0
  for v in iterator do
    t = t + 1
  end
  return t
end

function gff3_encode(s)
  return string.gsub(s, "[\t\n\r;=%&,]", function (c)
            return string.format("%%%02X", string.byte(c))
         end)
end

function gff3_decode(s)
  return string.gsub(s, "%%([0-9a-fA-F][1-9a-fA-F])", function (n)
            return string.char(tonumber("0x" .. n))
         end)
end

function gff3_extract_structure(str)
  local ret = {}
  for _,v in ipairs(str:split(",")) do
    local res = {}
    local v = gff3_decode(v)
    for _,pair in ipairs(v:split(";")) do
      key, value = unpack(pair:split("="))
      res[key] = value
    end
    table.insert(ret, res)
  end
  return ret
end

nodemt = debug.getregistry()["GenomeTools.genome_node"]
function nodemt.children_of_type(node, type)
  local nit = node:children()
  return function()
    local n = nit()
    while n and n:get_type() ~= type do
      n = nit()
    end
    return n
  end
end

function nodemt.children_of_supertype(node, type)
  local nit = node:children()
  return function()
    local n = nit()
    while n and not n:get_type():is_a(type) do
      n = nit()
    end
    return n
  end
end

function nodemt.children_matching_type(node, type_pat)
  local nit = node:children()
  return function()
    local n = nit()
    while n and string.match(n:get_type(), type_pat) do
      n = nit()
    end
    return n
  end
end

matchers = {
  should_be = function(value, expected)
    if value ~= expected then
      return false, "expecting "..tostring(expected)..", not ".. tostring(value)
    end
    return true
  end;

  should_be_truthy = function(value)
    if not value then
      return false, tostring(value) .. " is not truthy"
    end
    return true
  end;

  should_be_falsy = function(value)
    if value then
      return false, tostring(value) .. " is not falsy"
    end
    return true
  end;

  should_be_smaller_than = function(value, expected)
    if value >= expected then
      return false, tostring(value).." is larger than ".. tostring(expected)
    end
    return true
  end;

  should_be_larger_than = function(value, expected)
    if value <= expected then
      return false, tostring(value).." is smaller than ".. tostring(expected)
    end
    return true
  end;

  should_not_be = function(value, expected)
    if value == expected then
      return false, "should not be "..tostring(expected) .." but is ".. tostring(value)
    end
    return true
  end;

  should_have_key = function(value, expected)
    if value[expected] == nil then
      return false, tostring(value).." does not have key ".. tostring(expected)
    end
    return true
  end;

  should_not_have_key = function(value, expected)
    if value[expected] ~= nil then
      return false, tostring(value).." has key ".. tostring(expected)
    end
    return true
  end;

  should_error = function(f)
    if pcall(f) then
      return false, "expecting an error but received none"
    end
    return true
  end;

  should_match = function(value, pattern)
    if not string.find(value, pattern) then
      return false, value .. " does not match pattern "..pattern
    end
    return true
  end;

  should_not_match = function(value, pattern)
    if string.find(value, pattern) then
      return false, value .. " matches pattern "..pattern
    end
    return true
  end;

  should_contain = function(value, expected)
    if not table.contains(value, expected) then
      return false, table.pretty_array(value) .. " does not contain value ".. tostring(expected)
    end
    return true
  end;

  should_not_contain = function(value, expected)
    if table.contains(value, expected) then
      return false, table.pretty_array(value) .. " contains value ".. tostring(expected)
    end
    return true
  end;
}
matchers.should_equal = matchers.should_be

-- also make matchers available as expect(foo).to_X() instead of
-- expect(foo).should_X() to make them more similar to natural language
_matchers = {}
for m, f in pairs(matchers) do
  if m:match('^should_') then
    if m:match('^should_not_') then
      _matchers[m:gsub('should_not_', 'not_to_')] = f
    end
    _matchers[m:gsub('should_', 'to_')] = f
  end
  _matchers[m] = f
end
matchers = _matchers