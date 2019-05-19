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

lp = require 'cgilua/lp'

DocVisitorHTML = {}

local template_dir

function DocVisitorHTML:new(template_path, header)
  assert(template_path and header)
  template_dir = template_path
  o = {}
  setmetatable(o, self)
  self.__index = self
  o.header = header
  return o
end

local function include(template, env)
  assert(template)
  local template_path = template_dir .. template
  env = env or {}
  env.io = io
  env.os = os
  env.ipairs = ipairs
  return lp.include(template_path, env)
end

local function codify(str)
  assert(str)
  local res = string.gsub(str, "<<(.-)>>", "@@%1@@")
  res = string.gsub(res, "<(.-)>", "<code>%1</code>")
  res = string.gsub(res, "@@(.-)@@", "<code><%1></code>")
  res = string.gsub(res, " ([%a_][%a%d_%.]-%(%))", " <code>%1</code>")
  res = string.gsub(res, "___(.-)___", "<strong>%1</strong>")
  return string.gsub(res, "__(.-)__", "<em>%1</em>")
end

local function paragraphify(str)
  assert(str)
  return string.gsub(str, "\n\n", "</p><p>")
end

function DocVisitorHTML:show_header()
  include(self.header)
end

function DocVisitorHTML:visit_classes(classes)
  assert(classes)
  include("classes.lp", { classes = classes })
end

function DocVisitorHTML:visit_modules(modules)
  assert(modules)
  include("modules.lp", { modules = modules })
end

function DocVisitorHTML:visit_class(classname, comments)
  assert(classname)
  include("class.lp", { classname = classname })
  if comments then
    for i, _ in ipairs(comments) do
      comments[i] = paragraphify(codify(comments[i]))
    end
  include("class_comments.lp", { comments = comments })
  end
end

function DocVisitorHTML:visit_module(modulename)
  assert(modulename)
  include("module.lp", { modulename = modulename })
end

local sole_function_visited = false

function DocVisitorHTML:visit_sole_function(desc)
  if not sole_function_visited then
    include("sole_function.lp")
    sole_function_visited = true
  end
  self:visit_method(desc)
end

function DocVisitorHTML:visit_method(desc)
  assert(desc)
  local name
  local prototype = desc.name
  if desc.rval then
    name = desc.rval .. " " .. desc.name
  else
    name = desc.name
  end
  include("method.lp", { name = name, args = desc.args,
                         comment = codify(desc.comment),
                         prototype = prototype })
end

function DocVisitorHTML:visit_variable(desc)
  assert(desc)
  local name
  local prototype = desc.name
  if desc.type then
    name = desc.type .. " " .. desc.name
  else
    name = desc.name
  end
  include("variable.lp", { name = name,
                         comment = codify(desc.comment),
                         prototype = prototype })
end

function DocVisitorHTML:visit_funcdef(desc)
  assert(desc)
  include("funcdef.lp", { name = desc.name, comment = codify(desc.comment) })
end

function DocVisitorHTML:visit_index(names)
  assert(names)
  include("index.lp", { names = names })
end

function DocVisitorHTML:show_footer()
  include("footer.lp")
end
