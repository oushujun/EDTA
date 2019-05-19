#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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
require 'core/error'

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  extern "GtTypeChecker* gt_type_checker_ref(GtTypeChecker*)"
  extern "void gt_type_checker_is_valid_p(GtTypeChecker*, const char*, void*)"
  extern "void gt_type_checker_delete(GtTypeChecker*)"
  extern "GtTypeChecker* gt_type_checker_obo_new(const char*, GtError*)"
  extern "GtTypeChecker* gt_type_checker_builtin_new()"

  class TypeChecker
    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end

    def is_valid?(type)
      n = DL::malloc(GT::NATIVEINTSIZE)
      GT.gt_type_checker_is_valid_p(@checker, type, n)
      rval = (n[0, n.size].unpack("i")[0] == 1)
      return rval
    end

    def to_ptr
      @checker
    end
  end

  class TypeCheckerBuiltin < TypeChecker
    def initialize
      @checker = GT.gt_type_checker_builtin_new()
      @checker.free = GT::symbol("gt_type_checker_delete", "0P")
    end
  end

  class TypeCheckerOBO < TypeChecker
    def initialize(filename)
      if !File.exist?(filename) or File.directory?(filename) \
        or !File.readable?(filename)
        GT.gterror("invalid file: #{filename}")
      end
      err = GT::Error.new
      @checker = GT.gt_type_checker_obo_new(filename, err.to_ptr)
      if @checker == GT::NULL or @checker == nil then
        GT.gterror(err)
      else
        @checker.free = GT::symbol("gt_type_checker_delete", "0P")
      end
    end
  end
end
