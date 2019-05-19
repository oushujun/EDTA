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
from gt.extended.custom_stream import CustomStream
from gt.extended.custom_visitor_example import CustomVisitorExample
from gt.extended.feature_node import FeatureNode

# This is an example custom stream implementation in Python.
# It applies a CustomVisitorExample to each node read from the input stream.
# For example, it can be used like this:
#
#   ins = GFF3InStream("somefile.gff3")
#   cs = CustomStreamExample(ins)
#   outs = GFF3OutStream(cs)
#   outs.pull()
#
# Please raise a GTError if a user error in the stream is to be signaled.
# This will ensure that error handling is done properly in the underlying C
# code.


class CustomStreamExample(CustomStream):

    def __init__(self, instream):
        CustomStream.__init__(self)
        self.instream = instream
        self.visitor = CustomVisitorExample()

    def next(self):
        node = self.instream.next_tree()
        if node:
            node.accept(self.visitor)
        return node
