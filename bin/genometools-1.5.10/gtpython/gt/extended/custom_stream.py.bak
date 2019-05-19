#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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

from ctypes import CFUNCTYPE, c_void_p, c_int, POINTER
from gt.dlload import gtlib
from gt.core.error import Error, gterror, GTError
from gt.extended.genome_stream import GenomeStream

NextFunc = CFUNCTYPE(c_int, POINTER(c_void_p), c_void_p)
FreeFunc = CFUNCTYPE(c_void_p, c_void_p)


class CustomStream(GenomeStream):

    def __init__(self):

        try:
            self.next
        except AttributeError:
            gterror("%s does not implement mandatory method 'next'!"
                    % self.__class__.__name__)

        def next_w(nodepp, err):
            error = Error(err)
            try:
                nextnode = self.next()
                if nextnode:
                    nodepp[0] = nextnode.gn
                else:
                    nodepp[0] = None
                return 0
            except Error:
                import sys
                errmsg = sys.exc_info()[1]
                error.set(str(errmsg))
                return -1

        self.next_cb = NextFunc(next_w)

        def free_w(ptr):
            try:
                self.free()
            except AttributeError:
                pass

        self.free_cb = FreeFunc(free_w)
        self.gs = gtlib.gt_script_wrapper_stream_new(self.next_cb,
                                                     self.free_cb)
        self._as_parameter_ = self.gs

    def from_param(cls, obj):
        if not isinstance(obj, CustomStream):
            raise TypeError("argument must be a CustomStream")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, c_int, POINTER
        gtlib.gt_script_wrapper_stream_new.restype = c_void_p
        gtlib.gt_script_wrapper_stream_new.argtypes = [NextFunc, FreeFunc]

    register = classmethod(register)
