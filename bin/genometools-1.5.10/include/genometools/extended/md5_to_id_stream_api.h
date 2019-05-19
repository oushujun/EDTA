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

#ifndef MD5_TO_ID_STREAM_API_H
#define MD5_TO_ID_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

/* Implements the <GtNodeStream> interface. A <GtMD5ToIDStream> converts MD5
   fingerprints used as sequence IDs to ``regular'' ones. */
typedef struct GtMD5ToIDStream GtMD5ToIDStream;

/* Create a <GtMD5toIDStream*> which converts MD5 sequence IDs from nodes it
   retrieves from its <in_stream> to ``regular'' ones (with the help of the
   given <region_mapping>). Takes ownership of <region_mapping>! */
GtNodeStream* gt_md5_to_id_stream_new(GtNodeStream *in_stream,
                                      GtRegionMapping *region_mapping);

#endif
