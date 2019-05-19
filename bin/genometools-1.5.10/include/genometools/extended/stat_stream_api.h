/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STAT_STREAM_API_H
#define STAT_STREAM_API_H

#include <stdio.h>
#include "core/file_api.h"
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtStatStream> gathers statistics
   about the <GtGenomeNode> objects it retrieves from its node source and passes
   them along unmodified. */
typedef struct GtStatStream GtStatStream;

/* Create a <GtStatStream> object which gathers statistics about the
   <GtGenomeNode> objects it retrieves from its <in_stream> and returns them
   unmodified. Besides the basic statistics, statistics about the following
   distributions can be gathered, if the corresponding argument equals <true>:
   <gene_length_distribution>, <gene_score_distribution>,
   <exon_length_distribution>, <exon_number_distribution>,
   <intron_length_distribution>, <cds_length_distribution>.

   If <used_sources> equals <true>, it is recorded which source tags have been
   encountered. */
GtNodeStream* gt_stat_stream_new(GtNodeStream *in_stream,
                                 bool gene_length_distribution,
                                 bool gene_score_distribution,
                                 bool exon_length_distribution,
                                 bool exon_number_distribution,
                                 bool intron_length_distribution,
                                 bool cds_length_distribution,
                                 bool used_sources);
/* Write the statistics gathered by <stat_stream> to <outfp>. */
void          gt_stat_stream_show_stats(GtStatStream *stat_stream,
                                        GtFile *outfp);

#endif
