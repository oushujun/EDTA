/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTF_IN_STREAM_API_H
#define GTF_IN_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtGTFInStream> parses a GTF2.2
   file and returns it as a stream of <GtGenomeNode> objects. */
typedef struct GtGTFInStream GtGTFInStream;

/* Create a <GtGTFInStream*> which subsequently reads the GTF file with the
   given <filename>. If <filename> equals <NULL>, the GTF data is read from
   <stdin>. */
GtNodeStream* gt_gtf_in_stream_new(const char *filename);

#endif
