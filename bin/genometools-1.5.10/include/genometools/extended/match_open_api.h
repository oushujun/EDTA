/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_OPEN_API_H
#define MATCH_OPEN_API_H

 /* The <GtMatchOpen> class, implementing the <GtMatch> interface, is meant to
    store results in the OpenMatch format, e.g. as output by Vmatch. */
typedef struct GtMatchOpen GtMatchOpen;

#include "extended/match_api.h"

/* Creates a new <GtMatchOpen> object, storing long values <weight> in addition
   to the generic match contents <seqid1>, <seqid2>, <start_seq1>, <start_seq2>,
   <end_seq1>, and <end_seq2>. */
GtMatch* gt_match_open_new(char *seqid1,
                           char *seqid2,
                           GtUword start_seq1,
                           GtUword start_seq2,
                           GtUword end_seq1,
                           GtUword end_seq2,
                           GtWord weight,
                           GtMatchDirection dir);

/* Sets <weight> to be the weight value in <mo>. */
void gt_match_open_set_weight(GtMatchOpen *mo, GtWord weight);

/* Returns the weight value stored in <mo>. */
GtWord gt_match_open_get_weight(GtMatchOpen *mo);

#endif
