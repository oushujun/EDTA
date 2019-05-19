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

-- Show genome node on stdout (using the optional <gff3_visitor>).
function GenomeTools_genome_node:show(gff3_visitor)
  local gff3_visitor = gff3_visitor or gt.gff3_visitor_new()
  self:accept(gff3_visitor)
end

-- Show marked parts of genome node on stdout.
function GenomeTools_genome_node:show_marked()
  if self:contains_marked() then
    local gni = gt.genome_node_iterator_new(self)
    local gn = gni:next()
    while gn do
      if gn:is_marked() then
        gn:output_leading()
        print("")
      end
      gn = gni:next()
    end
  end
end
