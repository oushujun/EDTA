/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef LTR_CLASSIFY_STREAM_API_H
#define LTR_CLASSIFY_STREAM_API_H

#include "core/error_api.h"
#include "core/hashmap_api.h"
#include "extended/node_stream_api.h"

typedef struct GtLTRClassifyStream GtLTRClassifyStream;

GtNodeStream* gt_ltr_classify_stream_new(GtNodeStream *in_stream,
                                         GtHashmap *features,
                                         const char *famprefix,
                                         char **current_state,
                                         GtUword *progress,
                                         GtError *err);

#endif
