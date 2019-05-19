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

from .block import *
from .canvas import *
from .color import *
from .custom_track import *
from .diagram import *
from .feature_index import *
from .feature_stream import *
from .graphics import *
from .image_info import *
from .layout import *
from .rec_map import *
from .style import *

try:
    Block.register(gtlib)
    CanvasCairoFileBase.register(gtlib)
    CustomTrack.register(gtlib)
    Diagram.register(gtlib)
    # DiagramFromArray.register(gtlib)
    FeatureIndex.register(gtlib)
    Graphics.register(gtlib)
    GraphicsCairo.register(gtlib)
    ImageInfo.register(gtlib)
    Layout.register(gtlib)
    RecMap.register(gtlib)
    Style.register(gtlib)
except AttributeError:
    # fail gracefully when AnnotationSketch symbols are not present
    pass
