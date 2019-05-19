-- $Id: re.lua,v 1.32 2008/10/09 20:25:06 roberto Exp $

local m = require"lpeg"
local _G = _G
local tonumber, type, print, error = tonumber, type, print, error
local mt = getmetatable(m.P(0))

module "re"

local any = m.P(1)

-- Pre-defined names
Predef = { nl = m.P"\n" }


local I = m.P(function (s,i) print(i, s:sub(1, i-1)); return i end)


local function getdef (id, Defs)
  local c = Defs and Defs[id]
  if not c then error("undefined name: " .. id) end
  return c
end


local function patt_error (s, i)
  local msg = (#s < i + 20) and s:sub(i)
                             or s:sub(i,i+20) .. "..."
  msg = ("pattern error near '%s'"):format(msg)
  error(msg, 2)
end

local function mult (p, n)
  local np = m.P(true)
  while n >= 1 do
    if n%2 >= 1 then np = np * p end
    p = p * p
    n = n/2
  end
  return np
end

local function equalcap (s, i, c)
  if type(c) ~= "string" then return nil end
  local e = #c + i
  if s:sub(i, e - 1) == c then return e else return nil end
end


local S = (m.S(" \t\n") + "--" * (any - m.S"\n")^0)^0

local name = m.R("AZ", "az") * m.R("AZ", "az", "09")^0

local exp_follow = m.P"/" + ")" + "}" + ":}" + "~}" + name + -1

name = m.C(name)


-- identifiers only have meaning in a given environment
local Identifier = name * m.Carg(1)

local num = m.C(m.R"09"^1) * S / tonumber

local String = "'" * m.C((any - "'")^0) * "'" +
               '"' * m.C((any - '"')^0) * '"'


local Cat = "%" * Identifier / function (c,Defs)
  local cat =  Defs and Defs[c] or Predef[c]
  if not cat then error ("name '" .. c .. "' undefined") end
  return cat
end

local Range = m.Cs(any * (m.P"-"/"") * (any - "]")) / m.R

local item = Cat + Range + m.C(any)

local Class =
    "["
  * (m.C(m.P"^"^-1))    -- optional complement symbol
  * m.Cf(item * (item - "]")^0, mt.__add) /
                          function (c, p) return c == "^" and any - p or p end
  * "]"

local function adddef (t, k, Defs, exp)
  if t[k] then
    error("'"..k.."' already defined as a rule")
  else
    t[k] = exp
  end
  return t
end

local function firstdef (n, Defs, r) return adddef({n}, n, Defs, r) end



local exp = m.P{ "Exp",
  Exp = S * ( m.V"Grammar"
            + m.Cf(m.V"Seq" * ("/" * S * m.V"Seq")^0, mt.__add) );
  Seq = m.Cf(m.Cc(m.P"") * m.V"Prefix"^0 , mt.__mul)
        * (#exp_follow + patt_error);
  Prefix = "&" * S * m.V"Prefix" / mt.__len
         + "!" * S * m.V"Prefix" / mt.__unm
         + m.V"Suffix";
  Suffix = m.Cf(m.V"Primary" * S *
          ( ( m.P"+" * m.Cc(1, mt.__pow)
            + m.P"*" * m.Cc(0, mt.__pow)
            + m.P"?" * m.Cc(-1, mt.__pow)
            + "^" * ( m.Cg(num * m.Cc(mult))
                    + m.Cg(m.C(m.S"+-" * m.R"09"^1) * m.Cc(mt.__pow))
                    )
            + "->" * S * ( m.Cg(String * m.Cc(mt.__div))
                         + m.P"{}" * m.Cc(nil, m.Ct)
                         + m.Cg(Identifier / getdef * m.Cc(mt.__div))
                         )
            + "=>" * S * m.Cg(Identifier / getdef * m.Cc(m.Cmt))
            ) * S
          )^0, function (a,b,f) return f(a,b) end );
  Primary = "(" * m.V"Exp" * ")"
            + String / m.P
            + Class
            + Cat
            + "{:" * (name * ":" + m.Cc(nil)) * m.V"Exp" * ":}" /
                     function (n, p) return m.Cg(p, n) end
            + "=" * name / function (n) return m.Cmt(m.Cb(n), equalcap) end
            + m.P"{}" / m.Cp
            + "{~" * m.V"Exp" * "~}" / m.Cs
            + "{" * m.V"Exp" * "}" / m.C
            + m.P"." * m.Cc(any)
            + "<" * name * ">" / m.V;
  Definition = Identifier * S * '<-' * m.V"Exp";
  Grammar = m.Cf(m.V"Definition" / firstdef * m.Cg(m.V"Definition")^0, adddef) /
                m.P
}

local pattern = S * exp / m.P * (-any + patt_error)


function compile (p, defs)
  if m.type(p) == "pattern" then return p end   -- already compiled
  local cp = pattern:match(p, 1, defs)
  if not cp then error("incorrect pattern", 3) end
  return cp
end


local mem
local fmem
local gmem
local mt = {__mode = "v"}

function match (s, p, i)
  local cp = mem[p]
  if not cp then
    cp = compile(p)
    mem[p] = cp
  end
  return cp:match(s, i or 1)
end

function find (s, p, i)
  local cp = fmem[p]
  if not cp then
    cp = compile(p)
    cp = m.P{ m.Cp() * cp + 1 * m.V(1) }
    fmem[p] = cp
  end
  return cp:match(s, i or 1)
end

function gsub (s, p, rep)
  gmem[p] = gmem[p] or {}
  local cp = gmem[p][rep]
  if not cp then
    cp = compile(p)
    cp = m.Cs((cp / rep + 1)^0)
    gmem[p][rep] = cp
  end
  return cp:match(s)
end


function updatelocale ()
  m.locale(Predef)
  Predef.a = Predef.alpha
  Predef.c = Predef.cntrl
  Predef.d = Predef.digit
  Predef.g = Predef.graph
  Predef.l = Predef.lower
  Predef.p = Predef.punct
  Predef.s = Predef.space
  Predef.u = Predef.upper
  Predef.w = Predef.alnum
  Predef.x = Predef.xdigit
  Predef.A = any - Predef.a
  Predef.C = any - Predef.c
  Predef.D = any - Predef.d
  Predef.G = any - Predef.g
  Predef.L = any - Predef.l
  Predef.P = any - Predef.p
  Predef.S = any - Predef.s
  Predef.U = any - Predef.u
  Predef.W = any - Predef.w
  Predef.X = any - Predef.x
  mem = {}    -- restart memoization
  fmem = {}
  gmem = {}
  _G.setmetatable(mem, mt)
  _G.setmetatable(fmem, mt)
  _G.setmetatable(gmem, mt)
end


updatelocale()

