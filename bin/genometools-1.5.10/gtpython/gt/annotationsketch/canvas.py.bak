#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
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

from gt.dlload import gtlib
from gt.annotationsketch.style import Style
from gt.annotationsketch.image_info import ImageInfo
from gt.core.error import Error, gterror
from gt.core.gtstr import Str

GRAPHICS_PDF = 0
GRAPHICS_PNG = 1
GRAPHICS_PS = 2
GRAPHICS_SVG = 3


class Canvas:

    def __init__(self, *args):
        raise NotImplementedError("Please call the constructor of a Canvas " +
                                  "implementation.")

    def __del__(self):
        try:
            gtlib.gt_canvas_delete(self.canvas)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, Canvas):
            raise TypeError("argument must be a Canvas")
        return obj._as_parameter_

    from_param = classmethod(from_param)


class CanvasCairoFileBase(Canvas):

    def from_param(cls, obj):
        if not isinstance(obj, CanvasCairoFile):
            raise TypeError("argument must be a CanvasCairoFile")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def to_file(self, filename):
        err = Error()
        rval = gtlib.gt_canvas_cairo_file_to_file(self.canvas,
                                                  str(filename).encode(
                                                      'UTF-8'),
                                                  err)
        if rval != 0:
            gterror(err)

    def to_stream(self):
        from ctypes import string_at
        str = Str(None)
        gtlib.gt_canvas_cairo_file_to_stream(self.canvas, str._as_parameter_)
        return string_at(str.get_mem(), str.length())

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, c_ulong, c_int
        gtlib.gt_canvas_delete.restype = None
        gtlib.gt_canvas_delete.argtypes = [c_void_p]
        gtlib.gt_canvas_cairo_file_to_file.restype = c_int
        gtlib.gt_canvas_cairo_file_to_file.argtypes = [c_void_p,
                                                       c_char_p, c_void_p]
        gtlib.gt_canvas_cairo_file_to_stream.restype = c_char_p
        gtlib.gt_canvas_cairo_file_to_stream.argtypes = [c_void_p, c_void_p]
        gtlib.gt_canvas_cairo_file_new.restype = c_void_p
        gtlib.gt_canvas_cairo_file_new.argtypes = [c_void_p, c_int, c_ulong,
                                                   c_ulong, c_void_p, c_void_p]

    register = classmethod(register)


class CanvasCairoFile(CanvasCairoFileBase):

    def __init__(self, style, width, height, ii=None):
        Style.from_param(style)
        err = Error()
        iip = None
        if ii:
            iip = ii._as_parameter_
        canvas = gtlib.gt_canvas_cairo_file_new(style._as_parameter_,
                                                GRAPHICS_PNG, width, height, iip, err._as_parameter_)
        if canvas == None:
            gterror(err)
        self.canvas = canvas
        self._as_parameter_ = self.canvas


class CanvasCairoFilePNG(CanvasCairoFileBase):

    def __init__(self, style, width, height, ii=None):
        err = Error()
        iip = None
        if ii:
            iip = ii._as_parameter_
        canvas = gtlib.gt_canvas_cairo_file_new(style._as_parameter_,
                                                GRAPHICS_PNG, width, height, iip, err._as_parameter_)
        if canvas == None:
            gterror(err)
        self.canvas = canvas
        self._as_parameter_ = self.canvas


class CanvasCairoFilePDF(CanvasCairoFileBase):

    def __init__(self, style, width, height, ii=None):
        err = Error()
        iip = None
        if ii:
            iip = ii._as_parameter_
        canvas = gtlib.gt_canvas_cairo_file_new(style._as_parameter_,
                                                GRAPHICS_PDF, width, height, iip, err._as_parameter_)
        if canvas == None:
            gterror(err)
        self.canvas = canvas
        self._as_parameter_ = self.canvas


class CanvasCairoFilePS(CanvasCairoFileBase):

    def __init__(self, style, width, height, ii=None):
        err = Error()
        iip = None
        if ii:
            iip = ii._as_parameter_
        canvas = gtlib.gt_canvas_cairo_file_new(style._as_parameter_,
                                                GRAPHICS_PS, width, height, iip, err._as_parameter_)
        if canvas == None:
            gterror(err)
        self.canvas = canvas
        self._as_parameter_ = self.canvas


class CanvasCairoFileSVG(CanvasCairoFileBase):

    def __init__(self, style, width, height, ii=None):
        err = Error()
        iip = None
        if ii:
            iip = ii._as_parameter_
        canvas = gtlib.gt_canvas_cairo_file_new(style._as_parameter_,
                                                GRAPHICS_SVG, width, height, iip, err._as_parameter_)
        if canvas == None:
            gterror(err)
        self.canvas = canvas
        self._as_parameter_ = self.canvas
