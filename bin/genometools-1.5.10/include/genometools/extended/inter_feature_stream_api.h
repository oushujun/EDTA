/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef INTER_FEATURE_STREAM_API_H
#define INTER_FEATURE_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtInterFeatureStream> inserts new
   feature nodes between existing feature nodes of a certain type. */
typedef struct GtInterFeatureStream GtInterFeatureStream;

/* Create a <GtInterFeatureStream*> which inserts feature nodes of type
   <inter_type> between the feature nodes of type <outside_type> it retrieves
   from <in_stream> and returns them. */
GtNodeStream* gt_inter_feature_stream_new(GtNodeStream *in_stream,
                                          const char *outside_type,
                                          const char *inter_type);

#endif
