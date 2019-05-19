#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg
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
from gt.annotationsketch.canvas import Canvas
from gt.annotationsketch.diagram import Diagram
from gt.annotationsketch.style import Style
from gt.core.error import Error, gterror
from ctypes import c_ulong, c_void_p, c_int, c_char_p, CFUNCTYPE, byref, \
    POINTER, c_double, c_uint

TrackOrderingFunc = CFUNCTYPE(c_int, c_char_p, c_char_p, c_void_p)


class Layout:

    def __init__(self, diagram, width, style):
        err = Error()
        self.layout = gtlib.gt_layout_new(diagram._as_parameter_, width,
                                          style._as_parameter_, err._as_parameter_)
        if err.is_set():
            gterror(err)
        self._as_parameter_ = self.layout

    def __del__(self):
        try:
            gtlib.gt_layout_delete(self.layout)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, Layout):
            raise TypeError("argument must be a Layout")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def sketch(self, canvas):
        err = Error()
        had_err = gtlib.gt_layout_sketch(self.layout, canvas._as_parameter_,
                                         err._as_parameter_)
        if err.is_set():
            gterror(err)

    def get_height(self):
        err = Error()
        height = c_ulong()
        gtlib.gt_layout_get_height(self.layout, byref(height),
                                   err._as_parameter_)
        if err.is_set():
            gterror(err)
        return height.value

    def set_track_ordering_func(self, func):

        def trackorderer(s1, s2, data_ptr):
            ret = func(s1, s2)
            try:
                if ret is None:
                    raise ValueError
                else:
                    ret = int(ret)
            except ValueError:
                gterror("Track ordering function must return a number!")
            return ret

        self.tof_cb = TrackOrderingFunc(trackorderer)
        self.tof = trackorderer
        gtlib.gt_layout_set_track_ordering_func(self.layout, self.tof_cb)

    def register(cls, gtlib):
        gtlib.gt_layout_delete.restype = None
        gtlib.gt_layout_delete.argtypes = [c_void_p]
        gtlib.gt_layout_new.restype = c_void_p
        gtlib.gt_layout_new.argtypes = [c_void_p, c_uint, c_void_p, c_void_p]
        gtlib.gt_layout_sketch.restype = c_int
        gtlib.gt_layout_sketch.argtypes = [c_void_p, c_void_p, c_void_p]
        gtlib.gt_layout_set_track_ordering_func.argtypes = [c_void_p,
                                                            TrackOrderingFunc]
        gtlib.gt_layout_get_height.restype = c_int
        gtlib.gt_layout_get_height.argtypes = [c_void_p, POINTER(c_ulong),
                                               c_void_p]

    register = classmethod(register)
