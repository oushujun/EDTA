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
from gt.core.error import Error, gterror
from gt.extended.genome_node import GenomeNode
from ctypes import byref, c_void_p


class GenomeStream:

    def __init__(self, *args):
        raise NotImplementedError("Please call the constructor of a " +
                                  "GenomeStream implementation.")

    def __del__(self):
        try:
            gtlib.gt_node_stream_delete(self.gs)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, GenomeStream):
            raise TypeError("argument must be a GenomeStream")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def next_tree(self):
        err = Error()
        genome_node = c_void_p()
        rval = gtlib.gt_node_stream_next(self.gs, byref(genome_node),
                                         err._as_parameter_)
        if rval != 0:
            gterror(err)
        if genome_node.value == None:
            return None
        else:
            return GenomeNode.create_from_ptr(genome_node.value)

    def pull(self):
        err = Error()
        rval = gtlib.gt_node_stream_pull(self.gs, err._as_parameter_)
        if rval != 0:
            gterror(err)

    def register(cls, gtlib):
        from ctypes import c_int, c_void_p, POINTER
        gtlib.gt_node_stream_delete.argtypes = [c_void_p]
        gtlib.gt_node_stream_delete.restype = None
        gtlib.gt_node_stream_next.argtypes = [c_void_p, c_void_p, c_void_p]
        gtlib.gt_node_stream_next.restype = c_int
        gtlib.gt_node_stream_pull.argtypes = [c_void_p, c_void_p]
        gtlib.gt_node_stream_pull.restype = c_int

    register = classmethod(register)
