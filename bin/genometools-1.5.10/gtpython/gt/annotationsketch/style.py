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

from gt.dlload import gtlib, CollectFunc
from gt.annotationsketch.color import Color
from gt.core.error import Error, gterror
from gt.core.gtstr import Str
from gt.core.str_array import StrArray
from gt.extended.genome_node import GenomeNode

STYLE_OK = 0
STYLE_NOT_SET = 1
STYLE_ERROR = 2


class Style:

    def __init__(self, ptr=None):
        if ptr:
            self.style = ptr
            self.own = False
        else:
            e = Error()
            self.style = gtlib.gt_style_new(e._as_parameter_)
            if self.style == 0 or self.style == None:
                gterror(e)
            self.own = True
        self._as_parameter_ = self.style

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_style_delete(self.style)
            except AttributeError:
                pass

    def from_param(cls, obj):
        if not isinstance(obj, Style):
            raise TypeError("argument must be a Style")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def load_file(self, filename):
        err = Error()
        rval = gtlib.gt_style_load_file(self.style,
                                        str(filename).encode('UTF-8'), err)
        if rval != 0:
            gterror(err)

    def load_str(self, string):
        err = Error()
        strg = Str(str(string).encode("utf-8"))
        rval = gtlib.gt_style_load_str(self.style, strg, err)
        if rval != 0:
            gterror(err)

    def to_str(self):
        err = Error()
        string = Str()
        if gtlib.gt_style_to_str(self.style, string._as_parameter_,
                                 err._as_parameter_) == 0:
            return str(string)
        else:
            gterror(err)

    def clone(self):
        sty = Style()
        string = self.to_str()
        print(string)
        sty.load_str(str(string))
        return sty

    def get_color(self, section, key, gn=None):
        from ctypes import byref
        color = Color()
        err = Error()
        gnp = None
        if gn:
            gnp = gn._as_parameter_
        rval = gtlib.gt_style_get_color(self.style, section, key, byref(color),
                                        gnp, err._as_parameter_)
        if rval == STYLE_OK:
            return color
        elif rval == STYLE_NOT_SET:
            return None
        elif rval == STYLE_ERROR:
            gterror(err)

    def set_color(self, section, key, color):
        from ctypes import byref
        gtlib.gt_style_set_color(self.style, section, key, byref(color))

    def get_cstr(self, section, key, gn=None):
        string = Str()
        err = Error()
        gnp = None
        if gn:
            gnp = gn._as_parameter_
        rval = gtlib.gt_style_get_str(self.style, section, key,
                                      string._as_parameter_, gnp, err._as_parameter_)
        if rval == STYLE_OK:
            return str(string)
        elif rval == STYLE_NOT_SET:
            return None
        elif rval == STYLE_ERROR:
            gterror(err)

    def set_cstr(self, section, key, value):
        string = Str(str(value.encode("utf-8")))
        gtlib.gt_style_set_str(self.style, section, key, string._as_parameter_)

    def get_num(self, section, key, gn=None):
        from ctypes import c_double, byref
        double = c_double()
        err = Error()
        gnp = None
        if gn:
            gnp = gn._as_parameter_
        rval = gtlib.gt_style_get_num(self.style, section, key, byref(double),
                                      gnp, err._as_parameter_)
        if rval == STYLE_OK:
            return double.value
        elif rval == STYLE_NOT_SET:
            return None
        elif rval == STYLE_ERROR:
            gterror(err)

    def set_num(self, section, key, number):
        from ctypes import c_double
        num = c_double(number)
        gtlib.gt_style_set_num(self.style, section, key, num)

    def get_bool(self, section, key, gn=None):
        from ctypes import byref, c_int
        bool = c_int()
        err = Error()
        gnp = None
        if gn:
            gnp = gn._as_parameter_
        rval = gtlib.gt_style_get_bool(self.style, section, key, byref(bool),
                                       gnp, err._as_parameter_)
        if rval == STYLE_OK:
            if bool.value == 1:
                return True
            else:
                return False
        elif rval == STYLE_NOT_SET:
            return None
        elif rval == STYLE_ERROR:
            gterror(err)

    def set_bool(self, section, key, val):
        if val == True:
            gtlib.gt_style_set_bool(self.style, section, key, 1)
        else:
            gtlib.gt_style_set_bool(self.style, section, key, 0)

    def unset(self, section, key):
        gtlib.gt_style_unset(self.style, section, key)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_double, c_float, c_void_p, \
            POINTER, c_int
        gtlib.gt_style_delete.restype = None
        gtlib.gt_style_delete.argtypes = [c_void_p]
        gtlib.gt_style_get_bool.restype = c_int
        gtlib.gt_style_get_bool.argtypes = [c_void_p, c_char_p, c_char_p,
                                            POINTER(c_int), c_void_p, c_void_p]
        gtlib.gt_style_get_color.restype = c_int
        gtlib.gt_style_get_color.argtypes = [c_void_p, c_char_p,
                                             c_char_p, POINTER(Color), c_void_p, c_void_p]
        gtlib.gt_style_get_num.restype = c_int
        gtlib.gt_style_get_num.argtypes = [c_void_p, c_char_p, c_char_p,
                                           POINTER(c_double), c_void_p, c_void_p]
        gtlib.gt_style_get_str.restype = c_int
        gtlib.gt_style_get_str.argtypes = [c_void_p, c_char_p, c_char_p,
                                           c_void_p, c_void_p, c_void_p]
        gtlib.gt_style_load_file.restype = c_int
        gtlib.gt_style_load_file.argtypes = [c_void_p, c_char_p, c_void_p]
        gtlib.gt_style_load_str.restype = c_int
        gtlib.gt_style_load_str.argtypes = [c_void_p, c_void_p, c_void_p]
        gtlib.gt_style_new.restype = c_void_p
        gtlib.gt_style_new.argtypes = [c_void_p]
        gtlib.gt_style_set_bool.restype = None
        gtlib.gt_style_set_bool.argtypes = [
            c_void_p, c_char_p, c_char_p, c_int]
        gtlib.gt_style_set_color.restype = None
        gtlib.gt_style_set_color.argtypes = [c_void_p, c_char_p, c_char_p,
                                             POINTER(Color)]
        gtlib.gt_style_set_num.restype = None
        gtlib.gt_style_set_num.argtypes = [c_void_p, c_char_p, c_char_p,
                                           c_double]
        gtlib.gt_style_set_str.restype = None
        gtlib.gt_style_set_str.argtypes = [c_void_p, c_char_p, c_char_p,
                                           c_void_p]
        gtlib.gt_style_to_str.restype = int
        gtlib.gt_style_to_str.argtypes = [c_void_p, c_void_p, c_void_p]
        gtlib.gt_style_unset.restype = None
        gtlib.gt_style_unset.argtypes = [c_void_p, c_char_p, c_char_p]

    register = classmethod(register)
