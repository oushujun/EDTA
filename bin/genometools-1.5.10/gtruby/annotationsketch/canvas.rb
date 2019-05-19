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
require 'annotationsketch/graphics'
require 'core/str'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtCanvas* gt_canvas_cairo_file_new(GtStyle*, int, unsigned int, " + \
                                             "unsigned int, GtImageInfo*, " + \
                                             "GtError*)"
  extern "int gt_canvas_cairo_file_to_file(GtCanvasCairoFile*, const char*, " +
                                           "GtError*)"
  extern "int gt_canvas_cairo_file_to_stream(GtCanvasCairoFile*, GtStr*)"
  extern "void gt_canvas_delete(GtCanvas*)"
  
  class Canvas
    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end
  end

  class CanvasCairoFileBase < Canvas
    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end
  
    def to_file(filename)
      err = GT::Error.new()
      rval = GT.gt_canvas_cairo_file_to_file(@canvas, filename, err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def to_stream()
      str = GT::Str.new(nil)
      GT.gt_canvas_cairo_file_to_stream(@canvas, str.to_ptr)
      str.get_mem.to_s(str.length)
    end

    def to_ptr
      @canvas
    end
  end

  class CanvasCairoFile < CanvasCairoFileBase
    def initialize(style, width, height, ii = nil)
      err = GT::Error.new
      if ii.nil? then
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PNG, \
                                              width, height, GT::NULL, \
                                              err.to_ptr)
      else
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PNG, \
                                              width, height, ii.to_ptr, \
                                              err.to_ptr)
      end
      if @canvas == GT::NULL or @canvas == nil then
        GT.gterror(err)
      else
        @canvas.free = GT::symbol("gt_canvas_delete", "0P")
      end
    end
  end

  class CanvasCairoFilePNG < CanvasCairoFileBase
    def initialize(style, width, height, ii = nil)
      err = GT::Error.new
      if ii.nil? then
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PNG, \
                                              width, height, GT::NULL, \
                                              err.to_ptr)
      else
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PNG, \
                                              width, height, ii.to_ptr, \
                                              err.to_ptr)
      end
      if @canvas == GT::NULL or @canvas == nil then
        GT.gterror(err)
      else
        @canvas.free = GT::symbol("gt_canvas_delete", "0P")
      end
    end
  end

  class CanvasCairoFilePS < CanvasCairoFileBase
    def initialize(style, width, height, ii = nil)
      err = GT::Error.new
      if ii.nil? then
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PS, \
                                              width, height, GT::NULL, \
                                              err.to_ptr)
      else
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PS, \
                                              width, height, ii.to_ptr, \
                                              err.to_ptr)
      end
      if @canvas == GT::NULL or @canvas == nil then
        GT.gterror(err)
      else
        @canvas.free = GT::symbol("gt_canvas_delete", "0P")
      end
    end
  end

  class CanvasCairoFileSVG < CanvasCairoFileBase
    def initialize(style, width, height, ii = nil)
      err = GT::Error.new
      if ii.nil? then
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_SVG, \
                                              width, height, GT::NULL, \
                                              err.to_ptr)
      else
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_SVG, \
                                              width, height, ii.to_ptr, \
                                              err.to_ptr)
      end
      if @canvas == GT::NULL or @canvas == nil then
        GT.gterror(err)
      else
        @canvas.free = GT::symbol("gt_canvas_delete", "0P")
      end
    end
  end

  class CanvasCairoFilePDF < CanvasCairoFileBase
    def initialize(style, width, height, ii = nil)
     err = GT::Error.new
      if ii.nil? then
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PDF, \
                                              width, height, GT::NULL, \
                                              err.to_ptr)
      else
        @canvas = GT.gt_canvas_cairo_file_new(style.style, GT::GRAPHICS_PDF, \
                                              width, height, ii.to_ptr, \
                                              err.to_ptr)
      end
      if @canvas == GT::NULL or @canvas == nil then
        GT.gterror(err)
      else
        @canvas.free = GT::symbol("gt_canvas_delete", "0P")
      end
    end
  end
end
