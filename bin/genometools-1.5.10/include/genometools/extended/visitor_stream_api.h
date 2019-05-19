/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef VISITOR_STREAM_API_H
#define VISITOR_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/node_visitor_api.h"

/* Implements the <GtNodeStream> interface. */
typedef struct GtVisitorStream GtVisitorStream;

/* Create a new <GtVisitorStream*>, takes ownership of <node_visitor>.
   This stream applies <node_visitor> to each node which passes through it.
   Can be used to implement all streams with such a functionality. */
GtNodeStream* gt_visitor_stream_new(GtNodeStream *in_stream,
                                    GtNodeVisitor *node_visitor);

#endif
