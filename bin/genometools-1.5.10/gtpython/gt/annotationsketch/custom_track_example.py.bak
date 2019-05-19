#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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

from gt.annotationsketch.custom_track import CustomTrack
from gt.annotationsketch.color import Color
from gt.core.gtrange import Range


class CustomTrackExample(CustomTrack):

    def __init__(self):
        super(CustomTrackExample, self).__init__()

    def get_height(self):
        return 50

    def get_title(self):
        return "Sample track drawn by a Python script"

    def render(self, graphics, ypos, rng, style, error):
        from random import random
        data = []
        for i in range(0, 120):
            data.append(random())
        graphics.draw_curve_data(graphics.get_xmargins(), ypos, Color(0,
                                                                      0, 1, .6), data, 120, Range(0, 1), 40)
        return 0
