#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg
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

require 'gtdlload'
require 'gthelper'
require 'core/str'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "int gt_genome_node_accept(GtGenomeNode*, GenomeVisitor*, GtError*)"
  extern "GtGenomeNode* gt_genome_node_ref(GtGenomeNode*)"
  extern "GtStr*        gt_genome_node_get_seqid(GtGenomeNode*)"
  extern "unsigned long gt_genome_node_get_start(GtGenomeNode*)"
  extern "unsigned long gt_genome_node_get_end(GtGenomeNode*)"
  extern "const char* gt_genome_node_get_filename(GtGenomeNode*)"
  extern "unsigned int gt_genome_node_get_line_number(GtGenomeNode*)"
  extern "void gt_genome_node_delete(GtGenomeNode*)"

  class GenomeNode
    attr_reader :genome_node
    def initialize(node_ptr, newref=false)
      if newref then
        @genome_node = GT.gt_genome_node_ref(node_ptr)
      else
        @genome_node = node_ptr
      end
      @genome_node.free = GT::symbol("gt_genome_node_delete", "0P")
    end

    def ==(node)
      (node.genome_node == @genome_node)
    end

    def get_range
      (GT::gt_genome_node_get_start(@genome_node)..\
       GT::gt_genome_node_get_end(@genome_node))
    end

    def get_filename
      GT.gt_genome_node_get_filename(@genome_node)
    end

    def get_line_number
      GT.gt_genome_node_get_line_number(@genome_node)
    end

    def to_ptr
      @genome_node
    end

    def get_seqid
      str = GT::Str.new(GT.gt_genome_node_get_seqid(@genome_node))
      if !str.nil? then
        str.to_s
      else
        ""
      end
    end

    def accept(visitor)
      err = GT::Error.new()
      rval = GT.gt_genome_node_accept(@genome_node, visitor.to_ptr, err.to_ptr)
      if rval != 0
        GT.gterror(err)
      end
    end
  end
end
