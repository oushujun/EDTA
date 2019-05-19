/*
  Copyright (c) 2006-2007, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2007       Center for Bioinformatics, University of Hamburg

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

#ifndef SORT_STREAM_API_H
#define SORT_STREAM_API_H

#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtSortStream> sorts the
   <GtGenomeNode> objects it retrieves from its node source. */
typedef struct GtSortStream GtSortStream;

/* Create a <GtSortStream*> which sorts the genome nodes it retrieves from
   <in_stream> and returns them unmodified, but in sorted order. */
GtNodeStream* gt_sort_stream_new(GtNodeStream *in_stream);

#endif
