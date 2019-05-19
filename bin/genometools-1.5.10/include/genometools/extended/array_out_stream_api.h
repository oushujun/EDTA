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

#ifndef ARRAY_OUT_STREAM_API_H
#define ARRAY_OUT_STREAM_API_H

#include "core/array_api.h"
#include "core/error_api.h"
#include "extended/node_stream_api.h"

/* The <GtArrayOutStream> class implements the <GtNodeStream> interface.
   <GtArrayOutStream> takes nodes from an input stream and adds them to a
   <GtArray>. This stream can be used to obtain nodes for processing outside
   of the usual stream flow. */
typedef struct GtArrayOutStream GtArrayOutStream;

/* Creates a new <GtArrayInStream>, storing new references to <GtFeatureNode>s
   from <in_stream> in <nodes>. Note that the array must be set up to contain
   pointers to <GtGenomeNode>s! */
GtNodeStream* gt_array_out_stream_new(GtNodeStream *in_stream,
                                      GtArray *nodes, GtError *err);

/* Like <gt_array_out_stream_new()>, but not restricted to feature nodes. */
GtNodeStream* gt_array_out_stream_all_new(GtNodeStream *in_stream,
                                          GtArray *nodes, GtError *err);

#endif
