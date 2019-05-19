/*
  Copyright (c) 2007, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007       Center for Bioinformatics, University of Hamburg

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

#ifndef UNIQ_STREAM_API_H
#define UNIQ_STREAM_API_H

#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtUniqStream> filters out
   repeated features it retrieves from its node source. */
typedef struct GtUniqStream GtUniqStream;

/* Create a <GtUniqStream> object which filters out repeated feature node graphs
   it retrieves from the sorted <in_stream> and return all other nodes.  Two
   feature node graphs are considered to be __repeated__ if they have the same
   depth-first traversal and each corresponding feature node pair is similar
   according to the <gt_feature_node_is_similar()> method. For such a repeated
   feature node graph the one with the higher score (of the top-level feature)
   is kept. If only one of the feature node graphs has a defined score, this one
   is kept. */
GtNodeStream* gt_uniq_stream_new(GtNodeStream*);
#endif
