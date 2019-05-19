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
  extern "GtError* gt_error_new()"
  extern "const char* gt_error_get(const GtError*)"
  extern "ibool gt_error_is_set(const GtError*)"
  extern "const char* gt_error_get(const GtError*)"
  extern "void gt_error_set_nonvariadic(GtError*, const char*)"
  extern "void gt_error_unset(GtError*)"
  extern "void gt_error_delete(GtError*)"

  class Error
    def initialize(e = nil)
      if e.nil? then
        @error = GT.gt_error_new()
        @error.free = GT::symbol("gt_error_delete", "0P")
      else
        @error = e
      end
    end

    def get
      GT.gt_error_get(@error)
    end

    def set(str)
      GT.gt_error_set_nonvariadic(@error, str.to_s)
    end

    def is_set?
      GT.gt_error_is_set(@error)
    end

    def unset
      GT.gt_error_unset(@error)
    end

    def to_ptr
      @error
    end
  end
end
