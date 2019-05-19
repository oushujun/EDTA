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

from gt.dlload import gtlib
from gt.extended.genome_node import GenomeNode


class CommentNode(GenomeNode):

    def __init__(self):
        pass

    @classmethod
    def create_new(cls, comment):
        fn = gtlib.gt_comment_node_new(str(comment).encode("UTF-8"))
        n = cls.create_from_ptr(fn, True)
        return n

    def get_comment(self):
        return gtlib.gt_comment_node_get_comment(self.gn).decode("UTF-8")

    def register(cls, gtlib):
        from ctypes import c_char_p,  c_void_p
        gtlib.gt_comment_node_new.restype = c_void_p
        gtlib.gt_comment_node_new.argtypes = [c_char_p]
        gtlib.gt_comment_node_get_comment.restype = c_char_p
        gtlib.gt_comment_node_get_comment.argtypes = [c_void_p]

    register = classmethod(register)
