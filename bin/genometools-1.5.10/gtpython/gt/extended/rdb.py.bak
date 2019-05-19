#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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
from gt.core.error import Error, gterror


class RDB:

    def __init__(self, *args):
        raise NotImplementedError("Please call the constructor of a " +
                                  "RDB implementation.")

    def __del__(self):
        try:
            gtlib.gt_rdb_delete(self.rdb)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, RDB):
            raise TypeError("argument must be an RDB")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_void_p
        gtlib.gt_rdb_delete.restype = None
        gtlib.gt_rdb_delete.argtypes = [c_void_p]

    register = classmethod(register)


class RDBSqlite(RDB):

    def __init__(self, filename):
        err = Error()
        rdb = gtlib.gt_rdb_sqlite_new(filename, err._as_parameter_)
        if rdb == None:
            gterror(err)
        self.rdb = rdb
        self._as_parameter_ = self.rdb

    def from_param(cls, obj):
        if not isinstance(obj, RDBSqlite):
            raise TypeError("argument must be an RDBSqlite")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p
        gtlib.gt_rdb_sqlite_new.restype = c_void_p
        gtlib.gt_rdb_sqlite_new.argtypes = [c_char_p, c_void_p]

    register = classmethod(register)
