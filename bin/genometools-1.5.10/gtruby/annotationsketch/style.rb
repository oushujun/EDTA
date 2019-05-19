#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2008-2010 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg
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
require 'core/error'
require 'annotationsketch/color'

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  STYLE_OK      = 0
  STYLE_NOT_SET = 1
  STYLE_ERROR   = 2

  # looks weird, but apparently the way to pass proper pointers to external
  # functions
  DoubleArg = struct [
    "double val"
  ]

  typealias "bool", "ibool"

  # same as above
  BoolArg = struct [
    "bool val"
  ]

  extern "GtStyle* gt_style_new(GtError*)"
  extern "int gt_style_load_file(GtStyle*, const char*, GtError*)"
  extern "int gt_style_load_str(GtStyle*, GtStr*, GtError*)"
  extern "int gt_style_to_str(const GtStyle*, GtStr*, GtError*)"
  extern "int gt_style_get_color(GtStyle*, const char*, const char*, GtColor*, " +
                                 "GenomeNode*, GtError*)"
  extern "void gt_style_set_color(GtStyle*, const char*, const char*, GtColor*)"
  extern "int  gt_style_get_str(const GtStyle*, const char*, " +
                               "const char*, GtStr*, GtGenomeNode*, GtError*)"
  extern "void gt_style_set_str(GtStyle*, const char*, const char*, GtStr*)"
  extern "int  gt_style_get_num(const GtStyle*, const char*, " +
                               "const char*, double*, GtGenomeNode*, GtError*)"
  extern "void gt_style_set_num_p(GtStyle*, const char*, const char*, double*)"
  extern "int  gt_style_get_bool(const GtStyle*, const char*, " +
                                "const char*, bool*, GtGenomeNode*, GtError*)"
  extern "void gt_style_set_bool(GtStyle*, const char*, const char*, bool)"
  extern "void gt_style_unset(GtStyle*, const char*, const char*)"
  extern "void gt_style_delete(GtStyle*)"

  class Style
    attr_reader :style

    def initialize(s = nil)
      err = GT::Error.new()
      if s.nil? then
        @style = GT.gt_style_new(err.to_ptr)
        @style.free = GT::symbol("gt_style_delete", "0P")
      else
        @style = s
      end
      if not @style then GT.gterror(err) end

    end

    def load_file(filename)
      err = GT::Error.new()
      rval = GT.gt_style_load_file(@style, filename, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def load_str(str)
      err = GT::Error.new()
      str = GT::Str.new(str)
      rval = GT.gt_style_load_str(@style, str.to_ptr, err.to_ptr)
      if rval != 0 then
        GT.gterror(err)
      end
    end

    def to_str()
      err = GT::Error.new()
      str = GT::Str.new(nil)
      if GT.gt_style_to_str(@style, str.to_ptr, err.to_ptr) == 0
        str.to_s
      else
        GT.gterror(err)
      end
    end

    def clone
      sty = GT::Style.new()
      str = self.to_str()
      sty.load_str(str)
      sty
    end

    def get_color(section, key, gn = GT::NULL)
      color = GT::Color.malloc
      err = GT::Error.new
      case GT.gt_style_get_color(@style, section, key, color, gn, err.to_ptr)
        when STYLE_OK then
          color
        when STYLE_NOT_SET then
          nil
        when STYLE_ERROR then
          GT.gterror(err)
      end
    end

    def set_color(section, key, color)
      GT.gt_style_set_color(@style, section, key, color)
    end

    def get_cstr(section, key, gn = GT::NULL)
      str = GT::Str.new(nil)
      err = GT::Error.new
      case GT.gt_style_get_str(@style, section, key, str, gn, err.to_ptr)
        when STYLE_OK then
          str.to_s
        when STYLE_NOT_SET then
          nil
        when STYLE_ERROR then
          GT.gterror(err)
      end
    end

    def set_cstr(section, key, value)
      str = GT::Str.new(value)
      GT.gt_style_set_str(@style, section, key, str)
    end

    def get_num(section, key, gn = GT::NULL)
      double = DoubleArg.malloc
      err = GT::Error.new
      case GT.gt_style_get_num(@style, section, key, double, gn, err.to_ptr)
        when STYLE_OK then
          double.val
        when STYLE_NOT_SET then
          nil
        when STYLE_ERROR then
          GT.gterror(err)
      end
    end

    def set_num(section, key, number)
      num = number.to_f
      double = DoubleArg.malloc
      double.val = num
      GT.gt_style_set_num_p(@style, section, key, double)
    end

    def get_bool(section, key, gn = GT::NULL)
      bool = BoolArg.malloc
      err = GT::Error.new
      case GT.gt_style_get_bool(@style, section, key, bool, gn, err.to_ptr)
        when STYLE_OK then
          bool.val
        when STYLE_NOT_SET then
          nil
        when STYLE_ERROR then
          GT.gterror(err)
      end
    end

    def set_bool(section, key, val)
      GT.gt_style_set_bool(@style, section, key, val)
    end

    def unset(section, key)
      GT.gt_style_unset(@style, section, key)
    end

    def to_ptr
      @style
    end
  end
end
