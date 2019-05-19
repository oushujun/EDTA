/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef EXTRACT_FEATURE_STREAM_API_H
#define EXTRACT_FEATURE_STREAM_API_H

#include <stdio.h>
#include "core/file_api.h"
#include "core/trans_table_api.h"
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

/* Implements the <GtNodeStream> interface. A <GtExtractFeatureStream> extracts
   the corresponding sequences of features. */
typedef struct GtExtractFeatureStream GtExtractFeatureStream;

/* Create a <GtExtractFeatureStream*> which extracts the corresponding sequences
   of feature nodes (of the given <type>) it retrieves from <in_stream> and
   writes them in FASTA format (with the given <width>) to <outfp>. If <join> is
   <true>, features of the given <type> are joined together before the sequence
   is extracted. If <translate> is <true>, the sequences are translated into
   amino acid sequences before they are written to <outfp>. If <seqid> is <true>
   the sequence IDs of the extracted features are added to the FASTA header.
   If <target> is <true> the target IDs of the extracted features are added to
   the FASTA header. Takes ownership of <region_mapping>! */
GtNodeStream* gt_extract_feature_stream_new(GtNodeStream *in_stream,
                                            GtRegionMapping *region_Mapping,
                                            const char *type, bool join,
                                            bool translate, bool seqid,
                                            bool target, GtUword width,
                                            GtFile *outfp);

void          gt_extract_feature_stream_retain_id_attributes(
                                                       GtExtractFeatureStream*);

void          gt_extract_feature_stream_set_trans_table(
                                                       GtExtractFeatureStream*,
                                                       GtTransTable*);

void          gt_extract_feature_stream_show_coords(GtExtractFeatureStream*);

#endif
