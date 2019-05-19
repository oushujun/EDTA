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
from gt.core.error import gterror
from gt.core.gtrange import Range
from gt.core.gtstr import Str
from gt.extended.strand import strandchars
from gt.extended.feature_node import FeatureNode


class Block(object):

    def __init__(self, ptr):
        if ptr == 0 or ptr == None:
            gterror("Block pointer cannot be NULL (was: " + str(ptr) +
                    ")")
        self.block = gtlib.gt_block_ref(ptr)
        self._as_parameter_ = self.block

    def __del__(self):
        try:
            gtlib.gt_block_delete(self.block)
        except AttributeError:
            pass

    def get_range(self):
        r = gtlib.gt_block_get_range(self.block)
        return (r.start, r.end)

    def get_type(self):
        return gtlib.gt_block_get_type(self.block).decode('UTF-8')

    def has_only_one_fullsize_element(self):
        return gtlib.gt_block_has_only_one_fullsize_element(self.block) == \
            1

    def merge(self, block2):
        gtlib.gt_block_merge(self.block, block2._as_parameter_)

    def clone(self):
        return Block(gtlib.gt_block_clone(self.block))

    def get_caption(self):
        s = Str(gtlib.gt_block_get_caption(self.block))
        return str(s)

    def set_strand(self, strand):
        if not strand in strandchars:
            gterror("Invalid strand '%s' -- must be one of %s" % (strand,
                                                                  strandchars))
        gtlib.gt_block_set_strand(self.block, strandchars.index(strand))

    def get_strand(self):
        return strandchars[gtlib.gt_block_get_strand(self.block)]

    def get_top_level_feature(self):
        f = gtlib.gt_block_get_top_level_feature(self.block)
        if f != 0:
            return FeatureNode.create_from_ptr(f, True)
        else:
            return None

    def get_size(self):
        return int(gtlib.gt_block_get_size(self.block))

    def from_param(cls, obj):
        if not isinstance(obj, Block):
            raise TypeError("argument must be a Block")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, c_int, c_ulong
        gtlib.gt_block_caption_is_visible.restype = c_int
        gtlib.gt_block_caption_is_visible.argtypes = [c_void_p]
        gtlib.gt_block_clone.restype = c_void_p
        gtlib.gt_block_clone.argtypes = [c_void_p]
        gtlib.gt_block_delete.restype = None
        gtlib.gt_block_delete.argtypes = [c_void_p]
        gtlib.gt_block_get_caption.restype = c_void_p
        gtlib.gt_block_get_caption.argtypes = [c_void_p]
        gtlib.gt_block_get_range.restype = Range
        gtlib.gt_block_get_range.argtypes = [c_void_p]
        gtlib.gt_block_get_size.restype = c_ulong
        gtlib.gt_block_get_size.argtypes = [c_void_p]
        gtlib.gt_block_get_strand.restype = c_int
        gtlib.gt_block_get_strand.argtypes = [c_void_p]
        gtlib.gt_block_get_top_level_feature.restype = c_void_p
        gtlib.gt_block_get_top_level_feature.argtypes = [c_void_p]
        gtlib.gt_block_get_type.restype = c_char_p
        gtlib.gt_block_get_type.argtypes = [c_void_p]
        gtlib.gt_block_has_only_one_fullsize_element.restype = c_int
        gtlib.gt_block_has_only_one_fullsize_element.argtypes = [c_void_p]
        gtlib.gt_block_merge.restype = None
        gtlib.gt_block_merge.argtypes = [c_void_p, c_void_p]
        gtlib.gt_block_ref.restype = c_void_p
        gtlib.gt_block_ref.argtypes = [c_void_p]
        gtlib.gt_block_set_caption_visibility.restype = None
        gtlib.gt_block_set_caption_visibility.argtypes = [c_void_p,
                                                          c_int]
        gtlib.gt_block_set_strand.restype = None
        gtlib.gt_block_set_strand.argtypes = [c_void_p, c_int]

    register = classmethod(register)
