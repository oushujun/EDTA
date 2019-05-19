/*
  Copyright (c) 2006-2008, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008       Center for Bioinformatics, University of Hamburg

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

#ifndef CSA_STREAM_API_H
#define CSA_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"

#define GT_DEFAULT_JOIN_LENGTH 300

/* Implements the <GtNodeStream> interface. A <GtCSAStream> takes spliced
   alignments and transforms them into consensus spliced alignments. */
typedef struct GtCSAStream GtCSAStream;

/* Create a <GtCSAStream*> which takes spliced alignments from its <in_stream>
   (which are at most <join_length> many bases apart), transforms them into
   consensus spliced alignments, and returns them. */
GtNodeStream* gt_csa_stream_new(GtNodeStream *in_stream,
                                GtUword join_length);

#endif
