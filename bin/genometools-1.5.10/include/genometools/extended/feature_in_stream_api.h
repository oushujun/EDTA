/*
  Copyright (c) 2013-2014 Daniel S. Standage <daniel.standage@gmail.com>

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

#ifndef FEATURE_IN_STREAM_API_H
#define FEATURE_IN_STREAM_API_H

#include "extended/feature_index_api.h"
#include "extended/node_stream_api.h"

typedef struct GtFeatureInStream GtFeatureInStream;

/* Create a new <GtFeatureInStream> using the given <GtFeatureIndex> as the
   source of a node stream. */
GtNodeStream* gt_feature_in_stream_new(GtFeatureIndex *fi);

/* Instruct <stream> to deliver GtRegionNode objects whose ranges are
   specified such as with <##sequence-region> pragmas, rather than those
   inferred from the features. Specifically, use
   <gt_feature_index_get_orig_range_for_seqid> on the underlying feature
   index rather than the default <gt_feature_index_get_range_for_seqid>. */
void gt_feature_in_stream_use_orig_ranges(GtFeatureInStream *stream);

#endif
