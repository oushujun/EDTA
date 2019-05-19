--[[
  Copyright (c) 2008, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008       Center for Bioinformatics, University of Hamburg

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

require 'lpeg'
require 'fileutils'

DocParser = {}

-- Common Lexical Elements
local Any             = lpeg.P(1)
local Newline         = lpeg.P("\n")
local Whitespace      = lpeg.S(" \t\n")
local OptionalSpace   = Whitespace^0
local Space           = Whitespace^1
local Semicolon       = lpeg.P(";")
local DefineSeparator = lpeg.P("\\")

-- Lexical Elements of Lua
local LuaLongCommentStart  = lpeg.P("--[[")
local LuaLongCommentEnd    = lpeg.P("]]")
local LuaLongComment       = lpeg.Cc("long comment") *
                             lpeg.C(LuaLongCommentStart *
                                    (Any - LuaLongCommentEnd)^0 *
                                    LuaLongCommentEnd)
local LuaShortCommentStart = lpeg.P("--")
local LuaShortCommentEnd   = Newline
local LuaShortCommentLine  = OptionalSpace * LuaShortCommentStart *
                             lpeg.C((Any - LuaShortCommentEnd)^0) *
                             LuaShortCommentEnd
local LuaShortComment      = lpeg.Cc("comment") * LuaShortCommentLine^1
local LuaOptionalComment   = LuaShortComment +
                             lpeg.Cc("comment") * lpeg.Cc("undefined")
local LuaCommentStart      = LuaLongCommentStart + LuaShortCommentStart
local LuaEnd               = lpeg.P("end")
local LuaLocalFunction    =  lpeg.P("local") * OptionalSpace *
                             lpeg.P("function") * OptionalSpace *
                             lpeg.P(Any - lpeg.P("("))^1 *
                             lpeg.P("(") * (Any - lpeg.P(")"))^0 *
                             lpeg.P(")")
local LuaGlobalFunction    = -lpeg.P("local") * OptionalSpace *
                             lpeg.C(lpeg.P("function")) * OptionalSpace *
                             lpeg.Cc("") * lpeg.C(lpeg.P(Any - lpeg.P("("))^1) *
                             lpeg.P("(") * lpeg.C((Any - lpeg.P(")"))^0) *
                             lpeg.P(")")
local ExportLuaMethod      = lpeg.Ct(LuaOptionalComment * LuaGlobalFunction)
local CodeStop             = LuaCommentStart + LuaLocalFunction + LuaGlobalFunction
local LuaCode              = lpeg.Cc("code") * lpeg.C((Any - CodeStop)^1)

-- Lexical Elements of (Lua) C
local Character      = lpeg.R("AZ", "az") + lpeg.R("09") + lpeg.P("_") +
                       lpeg.P("*")
local CCommentStart  = lpeg.P("/*")
local CCommentEnd    = lpeg.P("*/")
local ExportLuaCComment = CCommentStart * lpeg.P(" exports the ") *
                          lpeg.Ct(lpeg.Cc("class") * lpeg.C(Character^1)) *
                          (Any - lpeg.P("to Lua:"))^1 *
                          lpeg.P("to Lua:") * ExportLuaMethod^0 *
                          (Any - CCommentEnd)^0 * CCommentEnd
local CComment       = CCommentStart * (Any - CCommentEnd)^0 * CCommentEnd
local CCode          = (Any - CCommentStart)^1

-- Lexical Elements of (pure) C
local Ifndef = lpeg.P("#ifndef") * Whitespace * Character^1 * Newline
local Define = lpeg.P("#define") * Whitespace * Character^1 * Newline
local Endif = lpeg.P("#endif") * Newline^0
local Include = lpeg.P("#include") * (Any - Newline)^1 * Newline
local ClassTypedef = lpeg.Ct(lpeg.Cc("class") *
                             (CCommentStart * lpeg.C((Any - CCommentEnd)^0) *
                               CCommentEnd)^0 * Newline^0 *
                             ((lpeg.P("typedef") * Space *
                               (lpeg.P("struct") + lpeg.P("enum")) * Space *
                               Character^1 * Space * lpeg.C(Character^1)) +
                              (lpeg.P("typedef") * Space * lpeg.P("char*") *
                               Space * lpeg.C(lpeg.P("GtTagValueMap")))) *
                              OptionalSpace * Semicolon)
local FunctionTypedef = lpeg.Ct(lpeg.Cc("funcdef") *
                                (CCommentStart * lpeg.C((Any - CCommentEnd)^0) *
                                CCommentEnd) * Newline^0 *
                                lpeg.P("typedef") * Space *
                                lpeg.C((Any - Semicolon)^1) * Semicolon)
local TypedefStruct = lpeg.P("typedef struct") * (Any - Semicolon)^1 * Semicolon
local OptionalWord = (Character^1 * Space)^-1
local Function = lpeg.Cc("function") *
                 lpeg.C(Character^1 * Space * OptionalWord * OptionalWord *
                        OptionalWord ) *
                 lpeg.C(lpeg.P(Any - lpeg.S("(;"))^1) * lpeg.P("(") *
                 lpeg.C((Any - lpeg.S(");"))^1) * lpeg.P(")") *
                 (Any - Semicolon)^0 * Semicolon
local FunctionPtr = lpeg.Cc("functionptr") *
                 lpeg.P("typedef") * Space *
                 lpeg.C(Character^1 * Space * OptionalWord * OptionalWord *
                        OptionalWord ) * OptionalSpace * lpeg.P("(") * lpeg.P("*") *
                 lpeg.C(lpeg.P(Any - lpeg.S("()"))^1) * lpeg.P(")") * lpeg.P("(") *
                 lpeg.C((Any - lpeg.P(")"))^1) * lpeg.P(")") *
                 (Any - Semicolon)^0 * Semicolon
local Variable = lpeg.Cc("variable") *
                 lpeg.C(lpeg.P("extern") * Space * Character^1 *Space*
                        OptionalWord * OptionalWord * OptionalWord) *
                 lpeg.C((Any - lpeg.S("();"))^0) * Semicolon
local ExportedComment = lpeg.Cc("comment") * CCommentStart *
                        lpeg.C((Any - CCommentEnd)^0) * CCommentEnd
local ExportedDefine = lpeg.Cc("function") *
                       (lpeg.P("#if") * (Any - Newline)^1 * Newline)^0 *
                       lpeg.C("#define") * Space *
                       lpeg.C(lpeg.P(Any - lpeg.P("("))^1) * lpeg.P("(") *
                       lpeg.C((Any - lpeg.P(")"))^1) * lpeg.P(")") *
                       OptionalSpace * DefineSeparator
local ExportedPlainDefine = lpeg.Cc("function") *
                            (lpeg.P("#if") * (Any - Newline)^1 * Newline)^0 *
                            lpeg.C("#define") * Space *
                            lpeg.C(lpeg.P(Any - (DefineSeparator + Space))^1) *
                            OptionalSpace * DefineSeparator
local ExportCMethod = lpeg.Ct(ExportedComment * Newline^0 * (Function + FunctionPtr + Variable))
local ExportCDefine = lpeg.Ct(ExportedComment * Newline^0 *
                              (ExportedDefine + ExportedPlainDefine))
local ModuleDef = lpeg.Ct(lpeg.Cc("module") * CCommentStart * Space *
                          lpeg.C(Character^1) * Space * lpeg.P("module") *
                          Space * CCommentEnd)

-- Lua Grammar
local Elem, Start = lpeg.V"Elem", lpeg.V"Start"
local LuaGrammar = lpeg.P{ Start,
  Start = lpeg.Ct(Elem^0);
  Elem  = ExportLuaMethod + LuaLongComment + LuaShortComment + Space +
          LuaLocalFunction + LuaCode;
}
LuaGrammar = LuaGrammar * -1

-- Lua C Grammar
local LuaCGrammar = lpeg.P{ Start,
 Start = lpeg.Ct(Elem^0);
 Elem  = lpeg.Ct(ExportLuaCComment) + CComment + Space + CCode;
}
LuaCGrammar = LuaCGrammar * -1

-- CGrammar
local CGrammar = lpeg.P{ Start,
  -- Start = lpeg.Ct(CComment * Newline^0 * Ifndef * Define * Elem^0 * Endif);
  Start = lpeg.Ct(CComment * Newline^0 * Ifndef * Define * Elem^0);
  Elem = ClassTypedef + ModuleDef + ExportCDefine + ExportCMethod + Space +
         Include + lpeg.C(TypedefStruct) + FunctionTypedef + CCode + CComment;
}
CGrammar = CGrammar * -1

function DocParser:new()
  o = {}
  o.lua_c_pattern = LuaCGrammar
  o.lua_pattern = LuaGrammar
  o.c_pattern = CGrammar
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocParser:parse(filename, be_verbose, is_lua)
  assert(filename)
  assert(is_header(filename) or is_lua_file(filename))
  if be_verbose then
    print("parsing " .. filename)
  end
  local file, err = io.open(filename, "r")
  assert(file, err)
  local filecontent = file:read("*a")
  if is_header(filename) then
    if is_lua then
      return lpeg.match(self.lua_c_pattern, filecontent)
    else
      return lpeg.match(self.c_pattern, filecontent)
    end
  else
    assert(is_lua_file(filename))
    return lpeg.match(self.lua_pattern, filecontent)
  end
end
