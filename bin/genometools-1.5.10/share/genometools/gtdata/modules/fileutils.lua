--[[
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

module(..., package.seeall)

require 'lfs'

-- returns true if file with <filename> exists, false otherwise
function file_exists(filename)
  assert(filename)
  if lfs.attributes(filename, "mode") then
    return true
  else
    return false
  end
end

local function has_mode(filename, mode)
  assert(filename and mode)
  local attr, err = lfs.attributes(filename, "mode")
  assert(attr, err)
  if attr == mode then
    return true
  else
    return false
  end
end

-- returns true if file with <filename> is a directory, false otherwise
function is_dir(filename)
  assert(filename)
  return has_mode(filename, "directory")
end

-- returns true if file with <filename> is a regular file, false otherwise
function is_regular_file(filename)
  assert(filename)
  return has_mode(filename, "file")
end

local function is_regular_file_with_ending(filename, ending)
  assert(filename and ending)
  local pattern = "%" .. ending .. "$"
  if string.find(filename, pattern) and is_regular_file(filename) then
    return true
  else
    return false
  end
end

-- returns true if file with <filename> is a header file, false otherwise
function is_header(filename)
  assert(filename)
  return is_regular_file_with_ending(filename, ".h")
end

-- returns true if file with <filename> is an API header file, false otherwise
function is_api_header(filename)
  assert(filename)
  return is_regular_file_with_ending(filename, "_api.h")
end

-- returns true if file with <filename> is a Lua file, false otherwise
function is_lua_file(filename)
  assert(filename)
  return is_regular_file_with_ending(filename, ".lua")
end
