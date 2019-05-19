#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "size_t", "unsigned long"
  extern "GtArray* gt_array_new(size_t)"
  extern "void* gt_array_get(const GtArray*, unsigned long)"
  extern "void  gt_array_add_ptr(GtArray*, void*)"
  extern "unsigned long gt_array_size(const GtArray*)"
  extern "void  gt_array_delete(GtArray*)"

  class Array
    def self.create(size = nil, own = true)
      if size.nil? then
        size = DL::sizeof("P")
      end
      return GT::Array.new(GT.gt_array_new(size), own)
    end

    def initialize(array_ptr, own = false)
      @array = array_ptr
      if own then
        @array.free = GT::symbol("gt_array_delete", "0P")
      end
    end

    def get(idx)
      GT.gt_array_get(@array, idx).ptr
    end

    def add(val)
      GT.gt_array_add_ptr(@array, val.to_ptr)
    end

    def size
      GT.gt_array_size(@array)
    end

    def to_ptr
      @array
    end
  end
end
