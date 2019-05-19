#
# Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

require 'gtdlload'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtStr* gt_str_new()"
  extern "GtStr* gt_str_new_cstr(const char*)"
  extern "void* gt_str_get_mem(const GtStr*)"
  extern "void gt_str_append_str(GtStr*, const GtStr*)"
  extern "void gt_str_append_cstr(GtStr*, const char*)"
  # we declare the return value as const char* instead of char*, because
  # otherwise dl/import wrongly assumes that it has responsibility for the
  # returned memory region (which leads to a double free())
  extern "const char* gt_str_get(const GtStr*)"
  extern "unsigned long gt_str_length(const GtStr*)"
  extern "void gt_str_delete(GtStr*)"

  class Str
    def initialize(cstr)
      if cstr.is_a?(String) then
        @str = GT.gt_str_new_cstr(cstr)
        @str.free = GT::symbol("gt_str_delete", "0P")
      elsif cstr.is_a?(DL::PtrData) then
        @str = cstr
      else
        @str = GT.gt_str_new()
        @str.free = GT::symbol("gt_str_delete", "0P")
      end

    end

    def append_str(str)
      GT.gt_str_append_str(@str, str)
    end

    def append_cstr(cstr)
      GT.gt_str_append_cstr(@str, cstr)
    end

    def to_s
      GT.gt_str_get(@str)
    end

    def get_mem
      GT.gt_str_get_mem(@str)
    end

    def to_ptr
      @str
    end

    def length
      GT.gt_str_length(@str)
    end
  end
end
