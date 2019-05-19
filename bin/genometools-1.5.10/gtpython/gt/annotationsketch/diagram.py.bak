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

from ctypes import CFUNCTYPE, c_char_p, c_void_p, c_int
from gt.dlload import gtlib
from gt.annotationsketch.block import Block
from gt.annotationsketch.canvas import Canvas
from gt.annotationsketch.custom_track import CustomTrack
from gt.annotationsketch.style import Style
from gt.core.array import Array
from gt.core.error import Error, gterror
from gt.core.gtrange import Range
from gt.core.gtstr import Str
from gt.extended.feature_index import FeatureIndex
from gt.extended.feature_node import FeatureNode

TrackSelectorFunc = CFUNCTYPE(c_void_p, c_void_p, c_void_p, c_void_p)


class Diagram:

    def from_array(arr, rng, style):
        from ctypes import byref, sizeof
        if rng.start > rng.end:
            gterror("range.start > range.end")
        gtarr = Array.create(sizeof(c_void_p), False)
        for i in arr:
            if not isinstance(i, FeatureNode):
                gterror("Diagram array must only contain FeatureNodes!")
            gtarr.add(i)
        diagram = gtlib.gt_diagram_new_from_array(gtarr._as_parameter_, byref(rng),
                                                  style._as_parameter_)
        return Diagram(diagram)

    from_array = staticmethod(from_array)

    def from_index(feature_index, seqid, rng, style):
        from ctypes import byref
        FeatureIndex.from_param(feature_index)
        Style.from_param(style)

        err = Error()
        if rng.start > rng.end:
            gterror("range.start > range.end")
        diagram = gtlib.gt_diagram_new(feature_index,
                                       str(seqid).encode('UTF-8'), byref(rng),
                                       style, err)
        if err.is_set():
            gterror(err)
        return Diagram(diagram)

    from_index = staticmethod(from_index)

    def __init__(self, ptr):
        self.diagram = ptr
        self._as_parameter_ = self.diagram

    def __del__(self):
        try:
            gtlib.gt_diagram_delete(self.diagram)
        except AttributeError:
            pass

    def set_track_selector_func(self, func):

        def trackselector(block_ptr, string_ptr, data_ptr):
            b = Block(block_ptr)
            string = Str(string_ptr)
            ret = func(b)
            if not ret:
                gterror("Track selector callback function must return a string!")
            string.append_cstr(ret)

        self.tsf_cb = TrackSelectorFunc(trackselector)
        self.tsf = trackselector
        gtlib.gt_diagram_set_track_selector_func(self.diagram, self.tsf_cb)

    def add_custom_track(self, ct):
        gtlib.gt_diagram_add_custom_track(self.diagram, ct._as_parameter_)

    def from_param(cls, obj):
        if not isinstance(obj, Diagram):
            raise TypeError("argument must be a Diagram")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, POINTER
        gtlib.gt_diagram_add_custom_track.restype = None
        gtlib.gt_diagram_add_custom_track.argtypes = [c_void_p,
                                                      c_void_p]
        gtlib.gt_diagram_delete.restype = None
        gtlib.gt_diagram_delete.argtypes = [c_void_p]
        gtlib.gt_diagram_new.restype = c_void_p
        gtlib.gt_diagram_new.argtypes = [c_void_p, c_char_p, POINTER(Range),
                                         c_void_p, c_void_p]
        gtlib.gt_diagram_new_from_array.restype = c_void_p
        gtlib.gt_diagram_new_from_array.argtypes = [c_void_p, POINTER(Range),
                                                    c_void_p]
        gtlib.gt_diagram_set_track_selector_func.restype = None
        gtlib.gt_diagram_set_track_selector_func.argtypes = [c_void_p,
                                                             TrackSelectorFunc]

    register = classmethod(register)
