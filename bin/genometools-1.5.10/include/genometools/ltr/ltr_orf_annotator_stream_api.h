/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef LTR_ORF_ANNOTATOR_STREAM_API_H
#define LTR_ORF_ANNOTATOR_STREAM_API_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "extended/node_stream_api.h"

/* implements the ``node stream'' interface */
typedef struct GtLTRORFAnnotatorStream GtLTRORFAnnotatorStream;

/* TODO: Should take a region mapping for easier ID->sequence mapping */
GtNodeStream* gt_ltr_orf_annotator_stream_new(GtNodeStream *in_stream,
                                              GtEncseq *encseq,
                                              unsigned int min,
                                              unsigned int max,
                                              bool all,
                                              GtError *err);

void          gt_ltr_orf_annotator_stream_set_progress_location(
                                                       GtLTRORFAnnotatorStream*,
                                                       GtUword*);
#endif
