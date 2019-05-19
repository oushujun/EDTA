/*
  Copyright (c) 2012 Sascha Kastens <mail@skastens.de>
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

#ifndef LTR_REFSEQ_MATCH_STREAM_API_H
#define LTR_REFSEQ_MATCH_STREAM_API_H

#include "core/error_api.h"
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

typedef struct GtLTRRefseqMatchStream GtLTRRefseqMatchStream;

/* Implements the <GtNodeStream> interface. <GtLTRRefseqMatchStream>
   extracts the sequences for each <GtGenomeNode> and matches them against
   a given set of reference sequences (<refseq_file>). If a match is found
   the associated <GtGenomeNode> will be extended by a new <GtFeatureNode> of
   type 'nucleotide_match' */
GtNodeStream* gt_ltr_refseq_match_stream_new(GtNodeStream *in_stream,
                                             const char *indexname,
                                             const char *refseq_file,
                                             const char *seq_file,
                                             double evalue,
                                             bool dust,
                                             int word_size,
                                             int gapopen,
                                             int gapextend,
                                             int penalty,
                                             int reward,
                                             int num_threads,
                                             double xdrop,
                                             double identity,
                                             const char *moreblast,
                                             bool flcands,
                                             double min_ali_len_perc,
                                             GtUword params_id,
                                             const char *source,
                                             GtError *err);

GtNodeStream* gt_ltr_refseq_match_stream_new_with_mapping(
                                             GtNodeStream *in_stream,
                                             const char *refseq_file,
                                             GtRegionMapping *rmap,
                                             double evalue,
                                             bool dust,
                                             int word_size,
                                             int gapopen,
                                             int gapextend,
                                             int penalty,
                                             int reward,
                                             int num_threads,
                                             double xdrop,
                                             double identity,
                                             const char *moreblast,
                                             bool flcands,
                                             double min_ali_len_perc,
                                             GtUword params_id,
                                             const char *source,
                                             GT_UNUSED GtError *err);

#endif
