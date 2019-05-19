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

#ifndef ARRAY_IN_STREAM_API_H
#define ARRAY_IN_STREAM_API_H

#include "core/array_api.h"
#include "core/error_api.h"
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. <GtArrayOutStream> takes
   an array of <GtGenomeNodes> and delivers them when used as an input stream.
   This stream can be used to feed nodes from outside into a stream flow. */
typedef struct GtArrayInStream GtArrayInStream;

/* Creates a new <GtArrayInStream>, delivering nodes from <nodes>. Note that
   the array must contain pointers to <GtGenomeNode>s! For every node passed,
   the value pointed to by <progress> is incremented by 1. */
GtNodeStream* gt_array_in_stream_new(GtArray *nodes, GtUword *progress,
                                     GtError *err);

#endif
