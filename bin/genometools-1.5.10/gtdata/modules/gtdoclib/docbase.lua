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

require 'stringext'
local w = require 'warning'

DocBase = {}

function DocBase:new()
  o = {}
  o.classes = {}
  o.classcomments = {}
  o.modules = {}
  o.moduledefs = {}
  o.variables = {}
  o.solefuncs = {}
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocBase:add_class(classname, comments, be_verbose)
  assert(classname)
  if be_verbose then
    print("class added: " .. classname)
  end
  self.classes[classname] = self.classes[classname] or {}
  self.classcomments[classname] = comments
end

function DocBase:add_module(modulename, be_verbose)
  assert(modulename)
  if be_verbose then
    print("module added: " .. modulename)
  end
  self.modules[modulename] = self.modules[modulename] or {}
  self.variables[modulename] = self.variables[modulename] or {}
end

function DocBase:add_method(funcret, funcname, funcargs, comment, be_verbose)
  assert(funcname and comment)
  local desc = {}
  -- remove ``GenomeTools_'' prefix which is used to extend exported C classes
  desc.rval = funcret
  desc.name = string.gsub(funcname, "^GenomeTools_", "")
  desc.args = funcargs
  desc.comment = comment
  if be_verbose then
    print("method added: " .. desc.name)
  end
  if self.last_module then
    self.modules[self.last_module][#self.modules[self.last_module] + 1] = desc
    return
  end
  local classname, match
  funcname = string.lower(string.gsub(desc.name, "_", ""))
  for class_to_search in pairs(self.classes) do
    local class_to_match = "^" .. string.lower(string.gsub(class_to_search, "_", ""))
    if be_verbose then
      print("match class: " .. class_to_match .. funcname)
    end
    match = string.match(funcname, class_to_match)
    if match then
      if not classname or string.len(match) > string.len(classname) then
        classname = class_to_search
      end
    end
  end
  if be_verbose and classname then
    print("classname found: " .. classname)
  end
  -- if this is a valid classname, try to store method in class
  if classname and self.classes[classname] then
    self.classes[classname][#self.classes[classname] + 1] = desc
  else
    self.solefuncs[#self.solefuncs + 1] = desc
  end
end

function DocBase:add_variable(vartype, varname, comment, be_verbose)
  assert(varname and comment)
  local desc = {}
  desc.type = vartype
  desc.name = varname
  desc.comment = comment
  if be_verbose then
    print("variable added: " .. desc.name)
  end
  if self.last_module then
    assert(self.variables)
    self.variables[self.last_module]
                  [#self.variables[self.last_module] + 1] = desc
    return
  end
  local classname, match
  varname = string.lower(string.gsub(desc.name, "_", ""))
  for class_to_search in pairs(self.classes) do
    local class_to_match = "^" .. string.lower(string.gsub(class_to_search, "_", ""))
    if be_verbose then
      print("match class: " .. class_to_match .. funcname)
    end
    match = string.match(funcname, class_to_match)
    if match then
      if not classname or string.len(match) > string.len(classname) then
        classname = class_to_search
      end
    end
  end
  if be_verbose and classname then
    print("classname found: " .. classname)
  end
  -- if this is a valid classname, try to store method in class
  if classname and self.classes[classname] then
    self.classes[classname][#self.classes[classname] + 1] = desc
  else
    self.solefuncs[#self.solefuncs + 1] = desc
  end
end

local function method_keyword(ast, be_verbose)
  for i, keyword in ipairs(ast) do
    if be_verbose then
      print("Try: " .. keyword)
    end
    if keyword == "function" or keyword == "functionptr" or keyword == "variable" then
      if be_verbose then
        print("Return: " .. i)
      end
      return i
    end
  end
  return 0
end

function DocBase:process_ast(ast, be_verbose)
  assert(ast)
  for _, v in ipairs(ast) do
    if type(v) == "table" then
      self:process_ast(v, be_verbose)
    else
      local keyword = ast[1]
      if be_verbose then
        print("keyword: " .. keyword)
      end
      if keyword == "class" then
        o.last_module = nil
        local comments
        if #ast > 2 then
          comments = {}
          for i = 2, #ast - 1 do
            if be_verbose then
              print("add class comment: " .. ast[i])
            end
            comments[#comments + 1] = ast[i]
          end
        end
        self["add_" .. ast[1]](self, ast[#ast], comments, be_verbose)
        break
      elseif keyword == "module" then
        self.last_module = ast[2]
        self["add_" .. ast[1]](self, ast[2], be_verbose)
      elseif keyword == "funcdef" then
        if be_verbose then
          print("funcdef keyword found")
        end
        if self.last_module then
          desc = {}
          desc.name = ast[3]
          desc.comment = ast[2]
          self.moduledefs[self.last_module] = self.moduledefs[self.last_module]
                                              or {}
          self.moduledefs[self.last_module][#self.moduledefs[self.last_module]
                                            + 1] = desc
        end
        break
      elseif keyword == "comment" then
        local funcpos = method_keyword(ast, be_verbose)
        local complete_comment = ""
        if funcpos > 0 then
          assert(funcpos > 2)
          assert(#ast == funcpos + 2 or #ast == funcpos + 3)
          if be_verbose then
            print("found: " .. ast[3] .. "!")
          end
          if ast[2] == "undefined" then
            w.warning("undefined comment")
          else
            complete_comment = table.concat(ast, "", 2, funcpos-1)
            complete_comment = string.strip(complete_comment)
          end
          if ast[3] == "variable" then
            self:add_variable(ast[funcpos+1], ast[funcpos+2], complete_comment,
                              be_verbose)
          else
            self:add_method(ast[funcpos+1], ast[funcpos+2], ast[funcpos+3],
                            complete_comment, be_verbose)
          end
          break
        elseif be_verbose then
          print("no function found!")
        end
      end
    end
  end
end

function DocBase:accept(visitor)
  assert(visitor)
  local method_names = {}
  -- visit all classes
  local sorted_classes = {}
  for classname in pairs(self.classes) do
    if #self.classes[classname] > 0 then
      sorted_classes[#sorted_classes + 1] = classname
    end
  end
  table.sort(sorted_classes)
  if visitor.visit_classes then
    visitor:visit_classes(sorted_classes)
  end
  -- visit all modules
  local sorted_modules = {}
  for modulename in pairs(self.modules) do
    if #self.modules[modulename] > 0 or
      #self.variables[modulename] > 0 or
      #self.moduledefs[modulename] > 0 then
      sorted_modules[#sorted_modules + 1] = modulename
    end
  end
  table.sort(sorted_modules)
  if visitor.visit_modules then
    visitor:visit_modules(sorted_modules)
  end
  -- visit sole functions
  for _, funcdesc in ipairs(self.solefuncs) do
    if visitor.visit_sole_function then
      visitor:visit_sole_function(funcdesc)
    else
      visitor:visit_method(funcdesc)
    end
    method_names[#method_names + 1] = funcdesc.name
  end
  -- visit each class
  for _, classname in ipairs(sorted_classes) do
    visitor:visit_class(classname, self.classcomments[classname])
    -- visit methods for class
    for _, method in ipairs(self.classes[classname]) do
      visitor:visit_method(method)
      method_names[#method_names + 1] = method.name
    end
  end
  -- visit each module
  if visitor.visit_module then
    for _, modulename in ipairs(sorted_modules) do
      visitor:visit_module(modulename)
      -- visit variables for module
      for _, variable in ipairs(self.variables[modulename]) do
        visitor:visit_variable(variable)
        method_names[#method_names + 1] = variable.name
      end
      -- visit funcdefs for module
      if self.moduledefs[modulename] then
        for _, funcdef in ipairs(self.moduledefs[modulename]) do
          visitor:visit_funcdef(funcdef)
        end
      end
      -- visit functions for module
      for _, method in ipairs(self.modules[modulename]) do
        visitor:visit_method(method)
        method_names[#method_names + 1] = method.name
      end
    end
  end
  -- visit all method and variable names (for index construction)
  if visitor.visit_index then
    table.sort(method_names)
    visitor:visit_index(method_names)
  end
end
