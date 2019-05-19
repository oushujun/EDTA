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

#ifndef CDS_STREAM_API_H
#define CDS_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

/* Implements the <GtNodeStream> interface. A <GtCDSStream> determines the
   coding sequence (CDS) for sequences determined by feature nodes of type
   __exon__ and adds them as feature nodes of type __CDS__. */
typedef struct GtCDSStream GtCDSStream;

/* Create a <GtCDSStream*> which determines the coding sequence (CDS) for
   sequences determined by feature nodes of type __exon__ it retrieves from
   <in_stream>, adds them as feature nodes of type __CDS__ and returns all
   nodes. <region_mapping> is used to map the sequence IDs of the feature nodes
   to the regions of the actual sequences. <minorflen> is the minimum length an
   ORF must have in order to be added. The CDS features are created with the
   given <source>. If <start_codon> equals <true> an ORF must begin with a start
   codon, otherwise it can start at any position. If <final_stop_codon> equals
   <true> the final ORF must end with a stop codon. If <generic_start_codons>
   equals <true>, the start codons of the standard translation scheme are used
   as start codons (otherwise the amino acid 'M' is regarded as a start codon).
*/
GtNodeStream* gt_cds_stream_new(GtNodeStream *in_stream,
                                GtRegionMapping *region_mapping,
                                unsigned int minorflen, const char *source,
                                bool start_codon, bool final_stop_codon,
                                bool generic_star_codons);

#endif
