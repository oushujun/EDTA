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
from gt.annotationsketch.rec_map import RecMap
import math


def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):

        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K


class ImageInfo:

    def __init__(self):
        self.ii = gtlib.gt_image_info_new()
        self._as_parameter_ = self.ii
        self.hotspots = None

    def __del__(self):
        try:
            gtlib.gt_image_info_delete(self.ii)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not (isinstance(obj, ImageInfo) or obj == None):
            raise TypeError("argument must be an ImageInfo")
        if obj == None:
            return None
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def get_height(self):
        return gtlib.gt_image_info_get_height(self.ii)

    def num_of_rec_maps(self):
        return gtlib.gt_image_info_num_of_rec_maps(self.ii)

    def compare_hotspots(cls, hs1, hs2):
        if hs1[2] - hs1[0] + 1 > hs2[2] - hs2[0] + 1:
            return 1
        elif hs1[2] - hs1[0] + 1 == hs2[2] - hs2[0] + 1:
            if hs1[3] > hs2[3]:
                return 1
            elif hs1[3] == hs2[3]:
                return 0
            else:
                return -1
        else:
            return -1

    compare_hotspots = classmethod(compare_hotspots)

    def each_hotspot(self):
        if not self.hotspots:
            self.hotspots = []
            for i in range(self.num_of_rec_maps()):
                rm = RecMap(gtlib.gt_image_info_get_rec_map(self.ii, i))
                self.hotspots.append([math.floor(rm.get_northwest_x()),
                                      math.floor(rm.get_northwest_y()),
                                      math.floor(rm.get_southeast_x()),
                                      math.floor(rm.get_southeast_y()),
                                      rm.get_genome_feature()])
            self.hotspots.sort(key=cmp_to_key(ImageInfo.compare_hotspots))
        for hs in self.hotspots:
            yield (hs[0], hs[1], hs[2], hs[3], hs[4])

    def register(cls, gtlib):
        from ctypes import c_void_p, c_ulong, c_uint
        gtlib.gt_image_info_delete.restype = None
        gtlib.gt_image_info_delete.argtypes = [c_void_p]
        gtlib.gt_image_info_get_rec_map.restype = c_void_p
        gtlib.gt_image_info_get_rec_map.argtypes = [c_void_p, c_ulong]
        gtlib.gt_image_info_num_of_rec_maps.restype = c_ulong
        gtlib.gt_image_info_num_of_rec_maps.argtypes = [c_void_p]
        gtlib.gt_image_info_get_height.restype = c_uint
        gtlib.gt_image_info_get_height.argtypes = [c_void_p]
        gtlib.gt_image_info_new.restype = c_void_p
        gtlib.gt_image_info_new.argtypes = []

    register = classmethod(register)
