#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

from gt.dlload import gtlib
from gt.extended.feature_node import FeatureNode


class RecMap:

    def __init__(self, rm):
        self.rm = rm
        self._as_parameter_ = self.rm

    def from_param(cls, obj):
        if not isinstance(obj, RecMap):
            raise TypeError("argument must be a RecMap")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def get_northwest_x(self):
        return gtlib.gt_rec_map_get_northwest_x(self.rm)

    def get_northwest_y(self):
        return gtlib.gt_rec_map_get_northwest_y(self.rm)

    def get_southeast_x(self):
        return gtlib.gt_rec_map_get_southeast_x(self.rm)

    def get_southeast_y(self):
        return gtlib.gt_rec_map_get_southeast_y(self.rm)

    def get_genome_feature(self):

        # refcount only this FeatureNode!

        return FeatureNode.create_from_ptr(gtlib.gt_rec_map_get_genome_feature(self.rm),
                                           True)

    def has_omitted_children(self):
        return gtlib.gt_rec_map_has_omitted_children(self.rm) == 1

    def register(cls, gtlib):
        from ctypes import c_int, c_void_p, c_double
        gtlib.gt_rec_map_get_northwest_x.restype = c_double
        gtlib.gt_rec_map_get_northwest_x.argtypes = [c_void_p]
        gtlib.gt_rec_map_get_northwest_y.restype = c_double
        gtlib.gt_rec_map_get_northwest_y.argtypes = [c_void_p]
        gtlib.gt_rec_map_get_southeast_x.restype = c_double
        gtlib.gt_rec_map_get_southeast_x.argtypes = [c_void_p]
        gtlib.gt_rec_map_get_southeast_y.restype = c_double
        gtlib.gt_rec_map_get_southeast_y.argtypes = [c_void_p]
        gtlib.gt_rec_map_get_genome_feature.restype = c_void_p
        gtlib.gt_rec_map_get_genome_feature.argtypes = [c_void_p]
        gtlib.gt_rec_map_has_omitted_children.restype = c_int
        gtlib.gt_rec_map_has_omitted_children.argtypes = [c_void_p]

    register = classmethod(register)
