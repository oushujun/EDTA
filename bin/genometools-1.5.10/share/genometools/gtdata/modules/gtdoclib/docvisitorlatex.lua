--[[
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

DocVisitorLaTeX = {}

local template_dir

function DocVisitorLaTeX:new(template_path, header)
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

local function trim(s)
  assert(s)
  return (string.gsub(s, "^%s*(.-)%s*$", "%1"))
end

local function align(s,count)
  local replace = string.rep(" ", count)
  return (replace .. string.gsub(s, ",%s+", ",\n" .. replace))
end

local function codify(str)
  if (str == nil) then
    str = ""
  end
  str = trim(str)
  local res = string.gsub(str, "\\", "\\textbackslash")
  res = string.gsub(res, "<<([a-zA-z_*()]-)>>", "@@%1@@")
  res = string.gsub(res, "<([a-zA-z_*()]-)>", "\\texttt{%1}")
  res = string.gsub(res, "@@([a-zA-z_*()]-)@@", "\\texttt{<%1>}")
  res = string.gsub(res, "%->", "\\(\\rightarrow\\)")
  res = string.gsub(res, "<%-", "\\(\\leftarrow\\)")
  res = string.gsub(res, ">=", "\\(\\ge\\)")
  res = string.gsub(res, ">", "\\(>\\)")
  res = string.gsub(res, "<=", "\\(\\le\\)")
  res = string.gsub(res, "<", "\\(<\\)")
  res = string.gsub(res, " ([%a_][%a%d_%.]-%(%))", " \\texttt{%1}")
  res = string.gsub(res, "___(.-)___", "\\textbf{%1}")
  res = string.gsub(res, "__(.-)__", "\\emph{%1}")
  res = string.gsub(res, "_", "\\_")
  res = string.gsub(res, "#", "\\#")
  res = string.gsub(res, "(%w+%^)%s*([%w\\{}]+)", "\\(%1{%2}\\)")
  res = string.gsub(res, "(%d+) ?\\%(", "\\(%1")
  res = string.gsub(res, "\\%) ?(%d+)", "%1\\)")
  res = string.gsub(res, "\\%) -%\\%(", " ")
  res = string.gsub(res, "%.%.%.", "\\ldots")
  return res
end

local function paragraphify(str)
  assert(str)
  return string.gsub(str, "\n\n", "\\\\")
end

function DocVisitorLaTeX:show_header()
  include(self.header)
end

function DocVisitorLaTeX:visit_classes(classes)
  assert(classes)
  include("classes_latex.lp", { classes = classes })
end

function DocVisitorLaTeX:visit_modules(modules)
  assert(modules)
  if (#modules == 0) then
    return
  end
  include("modules_latex.lp", { modules = modules })
end

function DocVisitorLaTeX:visit_class(classname, comments)
  assert(classname)
  include("class_latex.lp", { classname = codify(classname) })
  if comments then
    for i, _ in ipairs(comments) do
      comments[i] = paragraphify(codify(comments[i]))
    end
  include("class_comments_latex.lp", { comments = comments })
  end
end

function DocVisitorLaTeX:visit_module(modulename)
  assert(modulename)
  include("module_latex.lp", { modulename = modulename })
end

local sole_function_visited = false

function DocVisitorLaTeX:visit_sole_function(desc)
  if not sole_function_visited then
    include("sole_function_latex.lp")
    sole_function_visited = true
  end
  self:visit_method(desc)
end

function DocVisitorLaTeX:visit_method(desc)
  assert(desc)
  local name
  local args
  if desc.rval then
    name = trim(desc.rval) .. " " .. trim(desc.name)
  else
    name = trim(desc.name)
  end
  local aligncol = string.len(name) + 1
  if string.find(name, '_') then
    aligncol = string.find(name, '_') - 1
  end
  if desc.args then
    args = align(trim(desc.args),aligncol)
  else
    args=""
  end
  include("method_latex.lp",
          { name = name, args = args,
          comment = paragraphify(codify(desc.comment)) })
end

function DocVisitorLaTeX:visit_variable(desc)
  assert(desc)
  local name
  local prototype = desc.name
  if desc.type then
    name = desc.type .. " " .. desc.name
  else
    name = desc.name
  end
  include("variable_latex.lp", { name = codify(name),
                         comment = paragraphify(codify(desc.comment)) })
end

function DocVisitorLaTeX:visit_funcdef(desc)
  assert(desc)
  include("funcdef_latex.lp", { name = codify(desc.name),
                                comment = paragraphify(codify(desc.comment)) })
end

function DocVisitorLaTeX:visit_index(names)
end

function DocVisitorLaTeX:show_footer()
  include("footer_latex.lp")
end
