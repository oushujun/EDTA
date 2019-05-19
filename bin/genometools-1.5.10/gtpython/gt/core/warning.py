#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
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
from ctypes import CFUNCTYPE, c_void_p, c_char_p, addressof

funcdef = CFUNCTYPE(None, c_void_p, c_char_p, c_char_p)
defhand = funcdef.in_dll(gtlib, "gt_warning_default_handler")
gtlib.gt_warning_disable.restype = None
gtlib.gt_warning_disable.argtypes = []
gtlib.gt_warning_set_handler.restype = None
gtlib.gt_warning_set_handler.argtypes = [funcdef, c_void_p]


def warning_disable():
    gtlib.gt_warning_disable()


def warning_enable_default():
    gtlib.gt_warning_set_handler(addressof(defhand), None)
