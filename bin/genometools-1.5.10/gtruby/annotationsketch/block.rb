#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
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
require 'extended/strand'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtBlock*       gt_block_ref(GtBlock*)"
  extern "GtRange*       gt_block_get_range_ptr(const GtBlock*)"
  extern "const char*    gt_block_get_type(const GtBlock*)"
  extern "bool           gt_block_has_only_one_fullsize_element(const GtBlock*)"
  extern "void           gt_block_merge(GtBlock*, GtBlock*)"
  extern "GtBlock*       gt_block_clone(GtBlock*)"
  extern "void           gt_block_set_strand(GtBlock*, int)"
  extern "GtFeatureNode* gt_block_get_top_level_feature(const GtBlock*)"
  extern "int            gt_block_get_strand(const GtBlock*)"
  extern "unsigned long  gt_block_get_size(const GtBlock*)"
  extern "void           gt_block_delete(GtBlock*)"

  class Block
    def initialize(ptr)
      @block = GT.gt_block_ref(ptr)
      @block.free = GT::symbol("gt_block_delete", "0P")
    end

    def get_range
      ptr = GT.gt_block_get_range_ptr(@block)
      ptr.struct!("LL", :start, :stop)
      (ptr[:start]..ptr[:stop])
    end

    def has_only_one_fullsize_element
      GT.gt_block_has_only_one_fullsize_element(@block)
    end

    def get_type
      GT.gt_block_get_type(@block)
    end

    def merge(block)
      GT.gt_block_merge(@block, block)
    end

    def clone
      Block.new(GT.gt_block_clone(@block))
    end

    def set_strand(strand)
      if !GT::STRANDCHARS.include?(strand)
        GT::gterror("Invalid strand: '#{strand}'")
      end
      GT.gt_block_set_strand(@block, GT::STRANDCHARS.index(strand))
    end

    def get_strand
      GT::STRANDCHARS[GT.gt_block_get_strand(@block)]
    end

    def get_size
      GT::gt_block_get_size(@block)
    end

    def to_ptr
      @block
    end
  end
end
