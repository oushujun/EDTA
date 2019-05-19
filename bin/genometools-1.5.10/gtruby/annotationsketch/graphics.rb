#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy.to_f, modify.to_f, and distribute this software for any
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
require 'dl/struct'
require 'gthelper'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "void   gt_graphics_draw_text_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_text_clip_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_text_centered_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_text_right_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_colored_text_p(GtGraphics*, void*)"
  extern "double gt_graphics_get_text_height(GtGraphics*)"
  extern "double gt_graphics_get_text_width(GtGraphics*, const char*)"
  extern "double gt_graphics_get_image_width(GtGraphics*)"
  extern "double gt_graphics_get_image_height(GtGraphics*)"
  extern "void   gt_graphics_set_margins_p(GtGraphics*, void*)"
  extern "double gt_graphics_get_xmargins(GtGraphics*)"
  extern "double gt_graphics_get_ymargins(GtGraphics*)"
  extern "void   gt_graphics_draw_horizontal_line_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_vertical_line_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_line_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_box_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_dashes_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_caret_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_rectangle_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_arrowhead_p(GtGraphics*, void*)"
  extern "void   gt_graphics_draw_curve_data_p(GtGraphics*, void*, double*,
                                               unsigned long, unsigned long,
                                               unsigned long, unsigned long)"
  extern "int    gt_graphics_save_to_file(const GtGraphics*, const char*,
                                          GtError*)"
  extern "void   gt_graphics_save_to_stream(const GtGraphics*, GtStr*)"
  extern "void   gt_graphics_delete(GtGraphics*)"

  extern "GtGraphics*  gt_graphics_cairo_new(int, unsigned int, unsigned int)"

  GRAPHICS_PDF = 0
  GRAPHICS_PNG = 1
  GRAPHICS_PS  = 2
  GRAPHICS_SVG = 3

  ARROW_LEFT  = 0
  ARROW_RIGHT = 1
  ARROW_BOTH  = 2
  ARROW_NONE  = 3

  class Graphics
    def initialize(ptr)
      @g = ptr
    end

    def draw_text(x, y, text)
      ptr = [x.to_f, y.to_f, text].pack("DDp").to_ptr
      GT.gt_graphics_draw_text_p(@g, ptr)
    end

    def draw_text_centered(x, y, text)
      ptr = [x.to_f, y.to_f, text].pack("DDp").to_ptr
      GT.gt_graphics_draw_text_centered_p(@g, ptr)
    end

    def draw_text_right(x, y, text)
      ptr = [x.to_f, y.to_f, text].pack("DDp").to_ptr
      GT.gt_graphics_draw_text_right_p(@g, ptr)
    end

    def draw_colored_text(x, y, color, text)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, color[:r], color[:g], color[:b],               \
                color[:a], text]
      GT.gt_graphics_draw_colored_text_p(@g, params.pack("DDDDDDp").to_ptr)
    end

    def get_image_height
      GT.gt_graphics_get_image_height(@g)
    end

    def get_image_width
      GT.gt_graphics_get_image_width(@g)
    end

    def get_text_height
      GT.gt_graphics_get_text_height(@g)
    end

    def get_text_width(txt)
      GT.gt_graphics_get_text_width(@g, txt)
    end

    def set_margins(xmargs, ymargs)
      ptr = [xmargs.to_f, ymargs.to_f].pack("DD").to_ptr
      GT.gt_graphics_set_margins_p(@g, ptr)
    end

    def get_xmargins
      GT.gt_graphics_get_xmargins(@g)
    end

    def get_ymargins
      GT.gt_graphics_get_ymargins(@g)
    end

    def draw_line(x, y, xto, yto, color, stroke_width)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, xto.to_f, yto.to_f, color[:r], color[:g],      \
                color[:b], color[:a], stroke_width.to_f ]
      GT.gt_graphics_draw_line_p(@g, params.pack("DDDDDDDDD").to_ptr)
    end

    def draw_horizontal_line(x, y, color, width, stroke_width)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, color[:r], color[:g], color[:b], color[:a],    \
                width.to_f, stroke_width.to_f]
      GT.gt_graphics_draw_horizontal_line_p(@g, params.pack("DDDDDDDD").to_ptr)
    end

    def draw_vertical_line(x, y, color, len, stroke_width)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, color[:r], color[:g], color[:b], color[:a],    \
                len.to_f, stroke_width.to_f]
      GT.gt_graphics_draw_vertical_line_p(@g, params.pack("DDDDDDDD").to_ptr)
    end

    def draw_box(x, y, width, height, fillcolor, astatus, awidth, swidth,      \
                 scolor, dashed)
      fillcolor.struct!("DDDD", :r, :g, :b, :a)
      scolor.struct!("DDDD", :r, :g, :b, :a)
      if dashed then
        d_int = 1
      else
        d_int = 0
      end
      params = [x.to_f, y.to_f, width.to_f, height.to_f, fillcolor[:r],        \
                fillcolor[:g], fillcolor[:b], fillcolor[:a], awidth.to_f,      \
                swidth.to_f, scolor[:r],scolor[:g], scolor[:b], scolor[:a],    \
                d_int, astatus]
      GT.gt_graphics_draw_box_p(@g, params.pack("DDDDDDDDDDDDDDii").to_ptr)
    end

    def draw_dashes(x, y, width, height, astatus, awidth, swidth, scolor)
      scolor.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, width.to_f, height.to_f, awidth.to_f,          \
                swidth.to_f, scolor[:r],scolor[:g], scolor[:b], scolor[:a],    \
                astatus]
      GT.gt_graphics_draw_dashes_p(@g, params.pack("DDDDDDDDDDi").to_ptr)
    end

    def draw_caret(x, y, width, height, astatus, awidth, swidth, scolor)
      scolor.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, width.to_f, height.to_f, awidth, swidth,       \
                scolor[:r],scolor[:g], scolor[:b], scolor[:a], astatus]
      GT.gt_graphics_draw_caret_p(@g, params.pack("DDDDDDDDDDi").to_ptr)
    end

    def draw_rectangle(x, y, filled, fcolor, stroked, scolor, swidth,          \
                       width, height)
      scolor.struct!("DDDD", :r, :g, :b, :a)
      fcolor.struct!("DDDD", :r, :g, :b, :a)
      if filled then
        f_int = 1
      else
        f_int = 0
      end
      if stroked then
        s_int = 1
      else
        s_int = 0
      end
      params = [x.to_f, y.to_f, fcolor[:r], fcolor[:g], fcolor[:b], fcolor[:a],\
                scolor[:r], scolor[:g], scolor[:b], scolor[:a], swidth, width, \
                height, f_int.to_i, s_int.to_i ]
      GT.gt_graphics_draw_rectangle_p(@g, params.pack("DDDDDDDDDDDDDii").to_ptr)
    end

    def draw_arrowhead(x, y, color, status)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, color[:r], color[:g], color[:b], color[:a],    \
                status]
      GT.gt_graphics_draw_arrowhead_p(@g, params.pack("DDDDDDi").to_ptr)
    end

    def draw_curve_data(x, y, color, data, ndata, s, e, height)
      color.struct!("DDDD", :r, :g, :b, :a)
      params = [x.to_f, y.to_f, color[:r], color[:g], color[:b], color[:a]]
      GT.gt_graphics_draw_curve_data_p(@g, params.pack("DDDDDD").to_ptr,       \
                                       data.to_ptr, s, e, ndata, height)
    end

    def to_file(filename)
      err = GT::Error.new()
      if GT.gt_graphics_save_to_file(@g, filename, err) < 0 then
        GT::gterror(err)
      end
    end

    def to_stream
      str = GT::Str.new()
      GT.gt_graphics_save_to_stream(@g, str)
      str.get_mem.to_s(str.length)
    end

    def get_text_height
      GT.gt_graphics_get_text_height(@g)
    end

    def get_text_width(text)
      GT.gt_graphics_get_text_width(@g, text)
    end

    def to_ptr
      @g
    end
  end

  class GraphicsCairo < Graphics
  end

  class GraphicsCairoPDF < GraphicsCairo
    def initialize(width, height)
      @g = GT.gt_graphics_cairo_new(GT::GRAPHICS_PDF, width, height)
      @g.free = GT::symbol("gt_graphics_delete", "0P")
    end
  end

  class GraphicsCairoPNG < GraphicsCairo
    def initialize(width, height)
      @g = GT.gt_graphics_cairo_new(GT::GRAPHICS_PNG, width, height)
      @g.free = GT::symbol("gt_graphics_delete", "0P")
    end
  end

  class GraphicsCairoPS < GraphicsCairo
    def initialize(width, height)
      @g = GT.gt_graphics_cairo_new(GT::GRAPHICS_PS, width, height)
      @g.free = GT::symbol("gt_graphics_delete", "0P")
    end
  end

  class GraphicsCairoSVG < GraphicsCairo
    def initialize(width, height)
      @g = GT.gt_graphics_cairo_new(GT::GRAPHICS_SVG, width, height)
      @g.free = GT::symbol("gt_graphics_delete", "0P")
    end
  end
end
