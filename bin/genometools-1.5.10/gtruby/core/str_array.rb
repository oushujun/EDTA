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
  extern "GtStrArray* gt_str_array_new()"
  extern "void gt_str_array_add_cstr(GtStrArray*, const char*)"
  extern "const char* gt_str_array_get(const GtStrArray*, unsigned long)"
  extern "unsigned long gt_str_array_size(const GtStrArray*)"
  extern "void gt_str_array_delete(GtStrArray*)"

  class StrArray
    attr_reader :str_array
    def initialize(str_array_ptr = GT.gt_str_array_new(), own = true)
      @str_array = str_array_ptr
      if own != false then
        @str_array.free = GT::symbol("gt_str_array_delete", "0P")
      end
    end

    def add_list(list)
      list.each { |cstr| GT.gt_str_array_add_cstr(@str_array, cstr.to_s) }
    end

    def add(string)
      GT.gt_str_array_add_cstr(@str_array, string.to_s)
    end

    def to_a
      strings = []
      1.upto(GT.gt_str_array_size(@str_array)) do |i|
        strings.push(GT.gt_str_array_get(@str_array, i-1))
      end
      strings
    end

    def to_ptr
      @str_array
    end
  end
end
