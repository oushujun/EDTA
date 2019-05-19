/*
  Copyright (c) 2008, 2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef BED_IN_STREAM_API_H
#define BED_IN_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtBEDInStream> allows one to
   parse a BED file and return it as a stream of <GtGenomeNode> objects. */
typedef struct GtBEDInStream GtBEDInStream;

/* Return a <GtBEDInStream> object which subsequently reads the BED file with
   the given <filename>. If <filename> equals <NULL>, the BED data is read from
   <stdin>. */
GtNodeStream* gt_bed_in_stream_new(const char *filename);

/* Create BED features parsed by <bed_in_stream> with given <type> (instead of
   the default "BED_feature"). */
void          gt_bed_in_stream_set_feature_type(GtBEDInStream *bed_in_stream,
                                                const char *type);

/* Create thick BED features parsed by <bed_in_stream> with given <type>
   (instead of the default "BED_thick_feature"). */
void          gt_bed_in_stream_set_thick_feature_type(GtBEDInStream
                                                      *bed_in_stream,
                                                      const char *type);

/* Create BED blocks parsed by <bed_in_stream> with given <type> (instead of
   the default "BED_block"). */
void          gt_bed_in_stream_set_block_type(GtBEDInStream *bed_in_stream,
                                              const char *type);

#endif
