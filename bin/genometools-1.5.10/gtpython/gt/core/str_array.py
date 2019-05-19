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


class StrArray:

    def __init__(self, arr=None):
        if not arr:
            self.strarr = gtlib.gt_str_array_new()
        else:
            self.strarr = gtlib.gt_str_array_ref(arr)
        self._as_parameter_ = self.strarr

    def __del__(self):
        try:
            gtlib.gt_str_array_delete(self.strarr)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, StrArray):
            raise TypeError("argument must be a StrArray")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def to_list(self):
        result = []
        for i in range(gtlib.gt_str_array_size(self.strarr)):
            result.append(str(gtlib.gt_str_array_get(
                self.strarr, i).decode('UTF-8')))
        return result

    def size(self):
        return gtlib.gt_str_array_size(self.strarr)

    def get(self, i):
        return gtlib.gt_str_array_get(self.strarr, i).decode("UTF-8")

    def add(self, s):
        gtlib.gt_str_array_add_cstr(self.strarr, str(s).encode('UTF-8'))

    def register(cls, gtlib):
        from ctypes import c_void_p, c_char_p, c_ulong
        gtlib.gt_str_array_add_cstr.restype = None
        gtlib.gt_str_array_add_cstr.argtypes = [c_void_p, c_char_p]
        gtlib.gt_str_array_get.restype = c_char_p
        gtlib.gt_str_array_get.argtypes = [c_void_p, c_ulong]
        gtlib.gt_str_array_size.restype = c_ulong
        gtlib.gt_str_array_size.argtypes = [c_void_p]
        gtlib.gt_str_array_add_cstr.argtypes = [c_void_p, c_char_p]
        gtlib.gt_str_array_new.restype = c_void_p
        gtlib.gt_str_array_new.argtypes = []
        gtlib.gt_str_array_ref.restype = c_void_p
        gtlib.gt_str_array_ref.argtypes = [c_void_p]
        gtlib.gt_str_array_delete.restype = None
        gtlib.gt_str_array_delete.argtypes = [c_void_p]

    register = classmethod(register)
