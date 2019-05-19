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

local modname = ...
module(modname, package.seeall)

require 'gt'

-- all GenomeTools modules which should be loaded
local gtmodules = { "fileutils",
                    "gtlua.feature_index",
                    "gtlua.genome_features",
                    "gtlua.genome_node",
                    "gtlua.helper",
                    "gtlua.range" }

-- everything that will be exported to the gt table
local gtexport = {}

local function load_module(mod)
  assert(mod)
  local t = require(mod)
  for k, v in pairs(t) do
    if k ~= "_M" and k ~= "_NAME" and k~= "_PACKAGE" then
      assert(not gtexport[k]) -- symbol is undefined
      gtexport[k] = v -- record symbol for export
    end
  end
end

local function load_modules(modules)
  assert(modules)
  -- load all modules
  for _, mod in ipairs(modules) do
    load_module(mod)
  end
  -- export all symbols
  for k, v in pairs(gtexport) do
    assert(not gt[k]) -- symbol is undefined
    gt[k] = v -- export symbol
  end
end

-- Reload <gt> module.
function reload()
  -- remove all exported symbols from gt table
  for k in pairs(gtexport) do
    gt[k] = nil
  end
  -- mark all packages as unloaded
  for _, mod in ipairs(gtmodules) do
    package.loaded[mod] = nil
  end
  package.loaded[modname] = nil
  -- reload
  require(modname)
end

-- register reload() function in gt table
gt.reload = reload

load_modules(gtmodules)
