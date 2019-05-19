/*
  Copyright (c) 2009-2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef ADD_INTRONS_STREAM_API_H
#define ADD_INTRONS_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtAddIntronsStream> inserts new
   feature nodes with type __intron__ between existing feature nodes with type
   __exon__. This is a special case of the <GtInterFeatureStream>. */
typedef struct GtAddIntronsStream GtAddIntronsStream;

/* Create a <GtAddIntronsStream*> which inserts feature nodes of type __intron__
   between feature nodes of type __exon__ it retrieves from <in_stream> and
   returns them. */
GtNodeStream* gt_add_introns_stream_new(GtNodeStream *in_stream);

#endif
