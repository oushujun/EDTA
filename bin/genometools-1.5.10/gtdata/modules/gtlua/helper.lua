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

-- Export the content of <gt> table to the global environment.
function export()
  for k,v in pairs(gt) do
    _G[k] = v
  end
end

-- Call external 'display' program for file <filename>.
function display(filename)
  assert(filename and gt.file_exists(filename))
  if os.execute("display " .. filename) ~= 0 then
    io.stdout:write("\nexit (type 'y' to confirm)? ")
    if io.stdin:read() == "y" then
      print("bye")
      os.exit(0)
    end
  end
end

-- Show all keys and values of table <tbl>.
function show_table(tbl)
  assert(tbl)
  for k,v in pairs(tbl) do
    print(string.format("k=%s, v=%s", k, type(v)))
  end
end

-- Show content of the <gt> table.
function show(all)
  local a = {}
  for k in pairs(gt) do
    a[#a+1] = k
  end
  table.sort(a)
  for i,v in pairs(a) do
    if all then
      print(string.format("%s (%s)", v, type(gt[v])))
    else
      print(v)
    end
  end
end

-- Reload the <gt> module and export its content to the global environment.
function re()
  gt.reload()
  gt.export()
end
