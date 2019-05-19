#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

require 'dl/import'
require 'gthelper'
require "core/str"

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtRegionNode* gt_region_node_new(GtStr*, unsigned long,
                                            unsigned long)"

  class RegionNode < GenomeNode
    def self.create(seqid, start, stop)
      unless start <= stop
        raise(ArgumentError, "start (#{start}) > stop (#{stop})")
      else
        newrn = GT.gt_region_node_new(GT::Str.new(seqid.to_s), start, stop)
        return GT::RegionNode.new(newrn, true)
      end
    end

    def initialize(gn, newref=false)
      super(gn, newref)
    end
  end
end
