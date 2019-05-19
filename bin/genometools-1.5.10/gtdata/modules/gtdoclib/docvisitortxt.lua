--[[
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

DocVisitorTxt = {}

function DocVisitorTxt:new()
  o = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocVisitorTxt:visit_modules(modules)
  print("modules:")
  for _, mod in ipairs(modules) do
    print(mod)
  end
end

function DocVisitorTxt:visit_class(classname, comments)
  assert(classname)
  io.write(string.format("class: %s\n", classname))
  if comments then
    print("comments: " .. table.concat(comments))
  end
end

function DocVisitorTxt:visit_module(modulename)
  assert(modulename)
  io.write(string.format("module: %s\n", modulename))
end

function DocVisitorTxt:visit_method(desc)
  assert(desc)
  if desc.args then
    io.write(string.format("method:\n%s\n%s(%s)\n", desc.comment, desc.name,
             desc.args))
  else
    io.write(string.format("method:\n%s\n%s\n", desc.comment, desc.name))
           end
end

function DocVisitorTxt:visit_variable(desc)
  assert(desc)
  io.write(string.format("variable:\n%s\n%s %s\n", desc.comment, desc.type,
                         desc.name))
end

function DocVisitorTxt:visit_funcdef(desc)
  assert(desc)
  io.write(string.format("fundef:\n%s\n%s\n", desc.comment, desc.name))
end

function DocVisitorTxt:visit_index(names)
  assert(names)
  print("index:")
  for _, name in ipairs(names) do
    print(name)
  end
end
