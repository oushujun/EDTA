#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2008,2015 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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


class Str:

    def __init__(self, s=None):
        if s == None:
            self.strg = gtlib.gt_str_new()
            self.own = True
        elif isinstance(s, str):
            self.strg = gtlib.gt_str_new_cstr(str(s).encode('UTF-8'))
            self.own = True
        elif isinstance(s, bytes):
            self.strg = gtlib.gt_str_new_cstr(s)
            self.own = True
        else:
            self.strg = s
            self.own = False
        self._as_parameter_ = self.strg

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_str_delete(self.strg)
            except AttributeError:
                pass

    def __str__(self):
        return str(gtlib.gt_str_get(self.strg).decode('UTF-8'))

    def get(self):
        return gtlib.gt_str_get(self.strg).decode('UTF-8')

    def reset(self):
        gtlib.gt_str_reset(self.strg)

    def from_param(cls, obj):
        if not isinstance(obj, Str):
            raise TypeError("argument must be a Str")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def append_cstr(self, string):
        gtlib.gt_str_append_cstr(self.strg, str(string).encode('UTF-8'))

    def length(self):
        return gtlib.gt_str_length(self.strg)

    def get_mem(self):
        return gtlib.gt_str_get_mem(self.strg)

    def register(cls, gtlib):
        from ctypes import c_void_p, c_char_p, c_ulong
        gtlib.gt_str_new.restype = c_void_p
        gtlib.gt_str_new.argtypes = []
        gtlib.gt_str_new_cstr.restype = c_void_p
        gtlib.gt_str_new_cstr.argtypes = [c_char_p]
        gtlib.gt_str_append_cstr.restype = None
        gtlib.gt_str_append_cstr.argtypes = [c_void_p, c_char_p]
        gtlib.gt_str_get.restype = c_char_p
        gtlib.gt_str_get.argtypes = [c_void_p]
        gtlib.gt_str_get_mem.restype = c_void_p
        gtlib.gt_str_get_mem.argtypes = [c_void_p]
        gtlib.gt_str_length.restype = c_ulong
        gtlib.gt_str_length.argtypes = [c_void_p]
        gtlib.gt_str_reset.restype = None
        gtlib.gt_str_reset.argtypes = [c_void_p]
        gtlib.gt_str_delete.restype = None
        gtlib.gt_str_delete.argtypes = [c_void_p]

    register = classmethod(register)
