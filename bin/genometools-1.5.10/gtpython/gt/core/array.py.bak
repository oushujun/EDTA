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


class Array:

    def create(size=None, own=True):
        if size is None:
            from ctypes import c_void_p, sizeof
            size = sizeof(c_void_p)
        return Array(gtlib.gt_array_new(size), own)

    create = staticmethod(create)

    def __init__(self, arr, own=False):
        self.array = arr
        self._as_parameter_ = self.array
        self.own = own

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_array_delete(self.array)
            except AttributeError:
                pass

    def from_param(cls, obj):
        if not isinstance(obj, Array):
            raise TypeError("argument must be an Array")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def get(self, i):
        return gtlib.gt_array_get(self.array, i).contents

    def size(self):
        return gtlib.gt_array_size(self.array)

    def add(self, val):
        gtlib.gt_array_add_ptr(self.array, val._as_parameter_)

    def register(cls, gtlib):
        from ctypes import c_void_p, c_uint, c_ulong, POINTER
        gtlib.gt_array_new.restype = c_void_p
        gtlib.gt_array_new.argtypes = [c_uint]
        gtlib.gt_array_ref.restype = c_void_p
        gtlib.gt_array_ref.argtypes = [c_void_p]
        gtlib.gt_array_get.restype = POINTER(c_void_p)
        gtlib.gt_array_get.argtypes = [c_void_p, c_ulong]
        gtlib.gt_array_size.restype = c_ulong
        gtlib.gt_array_size.argtypes = [c_void_p]
        gtlib.gt_array_add_ptr.restype = None
        gtlib.gt_array_add_ptr.argtypes = [c_void_p, c_void_p]
        gtlib.gt_array_delete.restype = None
        gtlib.gt_array_delete.argtypes = [c_void_p]

    register = classmethod(register)
