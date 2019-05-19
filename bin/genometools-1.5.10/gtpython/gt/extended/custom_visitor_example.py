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

from gt.core.error import GTError
from gt.extended.custom_visitor import CustomVisitor
from gt.extended.feature_node import FeatureNode


# This is an example custom visitor implementation in Python.
# It adds a new subfeature to the visited node.
#
# For each node type, there is an accompanying visit_*_node method which is
# called when a node of this type is visited. If a method is missing, then the
# corresponding node type is ignored.
#
# Please raise a GTError if a user error in the visitor is to be signaled.
# This will ensure that error handling is done properly in the underlying C
# code.
class CustomVisitorExample(CustomVisitor):

    def __init__(self):
        CustomVisitor.__init__(self)

    def visit_feature_node(self, fn):
        new_child = FeatureNode.create_new(
            fn.get_seqid(), "bar", 100, 1000, "+")
        fn.add_child(new_child)

    def visit_region_node(self, rn):
        pass

    def visit_comment_node(self, cn):
        pass

    def visit_sequence_node(self, sn):
        pass
