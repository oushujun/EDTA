#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg
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

from gt.dlload import gtlib
from gt.annotationsketch.color import Color
from gt.core.array import Array
from gt.core.error import Error, gterror
from gt.core.gtrange import Range
from gt.core.gtstr import Str

GRAPHICS_PDF = 0
GRAPHICS_PNG = 1
GRAPHICS_PS = 2
GRAPHICS_SVG = 3

ARROW_LEFT = 0
ARROW_RIGHT = 1
ARROW_BOTH = 2
ARROW_NONE = 3


class Graphics:

    def __init__(self, p):
        self.g = p
        self.own = False
        self._as_parameter_ = self.g

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_graphics_delete(self.g)
            except AttributeError:
                pass

    def from_param(cls, obj):
        if not isinstance(obj, Graphics):
            raise TypeError("argument must be a Graphics")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def draw_text(self, x, y, text):
        gtlib.gt_graphics_draw_text(self.g, x, y, str(text).encode("utf-8"))

    def draw_text_centered(self, x, y, text):
        gtlib.gt_graphics_draw_text_centered(
            self.g, x, y, str(text).encode("utf-8"))

    def draw_text_right(self, x, y, text):
        gtlib.gt_graphics_draw_text_right(
            self.g, x, y, str(text).encode("utf-8"))

    def draw_colored_text(self, x, y, color, text):
        gtlib.gt_graphics_draw_colored_text(
            self.g, x, y, color, str(text).encode("utf-8"))

    def get_image_height(self):
        return gtlib.gt_graphics_get_image_height(self.g)

    def get_image_width(self):
        return gtlib.gt_graphics_get_image_width(self.g)

    def get_text_height(self):
        return gtlib.gt_graphics_get_text_height(self.g)

    def get_text_width(self, text):
        return gtlib.gt_graphics_get_text_width(self.g, str(text).encode("utf-8"))

    def set_margins(self, xmargs, ymargs):
        gtlib.gt_graphics_set_margins(self.g, xmargs, ymargs)

    def get_xmargins(self):
        return gtlib.gt_graphics_get_xmargins(self.g)

    def get_ymargins(self):
        return gtlib.gt_graphics_get_ymargins(self.g)

    def draw_line(self, x, y, xto, yto, color, stroke_width):
        gtlib.gt_graphics_draw_line(self.g, x, y, xto, yto, color,
                                    stroke_width)

    def draw_horizontal_line(self, x, y, color, width, stroke_width):
        gtlib.gt_graphics_draw_horizontal_line(self.g, x, y, color,
                                               width, stroke_width)

    def draw_vertical_line(self, x, y, color, length, stroke_width):
        gtlib.gt_graphics_draw_vertical_line(self.g, x, y, color, length,
                                             stroke_width)

    def draw_box(
        self,
        x,
        y,
        width,
        height,
        fillcolor,
        astatus,
        awidth,
        swidth,
        scolor,
        dashed,
    ):
        if dashed:
            d_int = 1
        else:
            d_int = 0
        gtlib.gt_graphics_draw_box(
            self.g,
            x,
            y,
            width,
            height,
            fillcolor,
            astatus,
            awidth,
            swidth,
            scolor,
            d_int,
        )

    def draw_dashes(
        self,
        x,
        y,
        width,
        height,
        astatus,
        awidth,
        swidth,
        scolor,
    ):
        gtlib.gt_graphics_draw_dashes(
            self.g,
            x,
            y,
            width,
            height,
            astatus,
            awidth,
            swidth,
            scolor,
        )

    def draw_caret(
        self,
        x,
        y,
        width,
        height,
        astatus,
        awidth,
        swidth,
        scolor,
    ):
        gtlib.gt_graphics_draw_caret(
            self.g,
            x,
            y,
            width,
            height,
            astatus,
            awidth,
            swidth,
            scolor,
        )

    def draw_rectangle(
        self,
        x,
        y,
        filled,
        fcolor,
        stroked,
        scolor,
        swidth,
        width,
        height,
    ):
        if filled:
            f_int = 1
        else:
            f_int = 0
        if stroked:
            s_int = 1
        else:
            s_int = 0
        gtlib.gt_graphics_draw_rectangle(
            self.g,
            x,
            y,
            f_int,
            fcolor,
            s_int,
            scolor,
            swidth,
            width,
            height,
        )

    def draw_arrowhead(self, x, y, color, status):
        gtlib.gt_graphics_draw_arrowhead(self.g, x, y, color, status)

    def draw_curve_data(self, x, y, color, data, ndata, valrange, height):
        from ctypes import c_double
        NDblArr = c_double * ndata
        cdata = NDblArr()
        for i in range(0, ndata):
            cdata[i] = data[i]
        gtlib.gt_graphics_draw_curve_data(self.g, x, y, color, cdata,
                                          ndata, valrange, height)

    def to_file(self, filename):
        err = Error()
        if gtlib.gt_graphics_save_to_file(self.g,
                                          str(filename).encode("utf-8"), err) < 0:
            gterror(err)

    def to_stream(self):
        from ctypes import string_at
        s = Str(None)
        gtlib.gt_graphics_save_to_stream(self.g, s._as_parameter_)
        return string_at(s.get_mem(), s.length())

    def get_text_height(self):
        return gtlib.gt_graphics_get_text_height(self.g)

    def get_text_width(self, text):
        return gtlib.gt_graphics_get_text_width(self.g, str(text).encode("utf-8"))

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, c_int, POINTER, c_double, \
            c_ulong
        gtlib.gt_graphics_draw_text.restype = None
        gtlib.gt_graphics_draw_text.argtypes = [c_void_p, c_double,
                                                c_double, c_char_p]
        gtlib.gt_graphics_draw_text_right.argtypes = [c_void_p, c_double,
                                                      c_double, c_char_p]
        gtlib.gt_graphics_draw_text_centered.argtypes = [c_void_p,
                                                         c_double, c_double, c_char_p]
        gtlib.gt_graphics_draw_colored_text.argtypes = [c_void_p,
                                                        c_double, c_double, Color, c_char_p]
        gtlib.gt_graphics_get_image_height.restype = c_double
        gtlib.gt_graphics_get_image_height.argtypes = [c_void_p]
        gtlib.gt_graphics_get_image_width.restype = c_double
        gtlib.gt_graphics_get_image_width.argtypes = [c_void_p]
        gtlib.gt_graphics_get_text_height.restype = c_double
        gtlib.gt_graphics_get_text_height.argtypes = [c_void_p]
        gtlib.gt_graphics_get_text_width.restype = c_double
        gtlib.gt_graphics_get_text_width.argtypes = [c_void_p, c_char_p]
        gtlib.gt_graphics_set_margins.argtypes = [c_void_p, c_double,
                                                  c_double]
        gtlib.gt_graphics_get_xmargins.restype = c_double
        gtlib.gt_graphics_get_xmargins.argtypes = [c_void_p]
        gtlib.gt_graphics_get_ymargins.restype = c_double
        gtlib.gt_graphics_get_ymargins.argtypes = [c_void_p]
        gtlib.gt_graphics_draw_vertical_line.argtypes = [c_void_p,
                                                         c_double, c_double, Color, c_double, c_double]
        gtlib.gt_graphics_draw_line.argtypes = [c_void_p, c_double,
                                                c_double, c_double, c_double, Color, c_double]
        gtlib.gt_graphics_draw_horizontal_line.argtypes = [c_void_p,
                                                           c_double, c_double, Color, c_double, c_double]
        gtlib.gt_graphics_draw_box.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_double,
            c_double,
            Color,
            c_int,
            c_double,
            c_double,
            Color,
            c_int,
        ]
        gtlib.gt_graphics_draw_dashes.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_double,
            c_double,
            c_int,
            c_double,
            c_double,
            Color,
        ]
        gtlib.gt_graphics_draw_caret.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_double,
            c_double,
            c_int,
            c_double,
            c_double,
            Color,
        ]
        gtlib.gt_graphics_draw_rectangle.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_int,
            Color,
            c_int,
            Color,
            c_double,
            c_double,
            c_double,
        ]
        gtlib.gt_graphics_draw_arrowhead.argtypes = [c_void_p, c_double,
                                                     c_double, Color, c_int]
        gtlib.gt_graphics_draw_curve_data.argtypes = [c_void_p, c_double,
                                                      c_double, Color, c_void_p, c_ulong, Range, c_ulong]
        gtlib.gt_graphics_save_to_file.restype = c_int
        gtlib.gt_graphics_save_to_file.argtypes = [c_void_p, c_char_p,
                                                   c_void_p]
        gtlib.gt_graphics_save_to_stream.argtypes = [c_void_p, c_void_p]
        gtlib.gt_graphics_delete.restype = None
        gtlib.gt_graphics_delete.argtypes = [c_void_p]

    register = classmethod(register)


class GraphicsCairo(Graphics):

    def __init__(self, p):
        self.g = p
        self.own = False
        self._as_parameter_ = self.g

    def from_param(cls, obj):
        if not isinstance(obj, GraphicsCairo):
            raise TypeError("argument must be a GraphicsCairo")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_ulong, c_void_p, c_int
        gtlib.gt_graphics_cairo_new.restype = c_void_p
        gtlib.gt_graphics_cairo_new.argtypes = [c_int, c_ulong, c_ulong]

    register = classmethod(register)


class GraphicsCairoPNG(GraphicsCairo):

    def __init__(self, width, height):
        self.g = gtlib.gt_graphics_cairo_new(GRAPHICS_PNG, width, height)
        self.own = True
        self._as_parameter_ = self.g

    def from_param(cls, obj):
        if not isinstance(obj, GraphicsCairoPNG):
            raise TypeError("argument must be a GraphicsCairoPNG")
        return obj._as_parameter_

    from_param = classmethod(from_param)


class GraphicsCairoPDF(GraphicsCairo):

    def __init__(self, width, height):
        self.g = gtlib.gt_graphics_cairo_new(GRAPHICS_PDF, width, height)
        self.own = True
        self._as_parameter_ = self.g

    def from_param(cls, obj):
        if not isinstance(obj, GraphicsCairoPDF):
            raise TypeError("argument must be a GraphicsCairoPDF")
        return obj._as_parameter_

    from_param = classmethod(from_param)


class GraphicsCairoPS(GraphicsCairo):

    def __init__(self, width, height):
        self.g = gtlib.gt_graphics_cairo_new(GRAPHICS_PS, width, height)
        self.own = True
        self._as_parameter_ = self.g

    def from_param(cls, obj):
        if not isinstance(obj, GraphicsCairoPS):
            raise TypeError("argument must be a GraphicsCairoPS")
        return obj._as_parameter_

    from_param = classmethod(from_param)


class GraphicsCairoSVG(GraphicsCairo):

    def __init__(self, width, height):
        self.g = gtlib.gt_graphics_cairo_new(GRAPHICS_SVG, width, height)
        self.own = True
        self._as_parameter_ = self.g

    def from_param(cls, obj):
        if not isinstance(obj, GraphicsCairoSVG):
            raise TypeError("argument must be a GraphicsCairoSVG")
        return obj._as_parameter_

    from_param = classmethod(from_param)
