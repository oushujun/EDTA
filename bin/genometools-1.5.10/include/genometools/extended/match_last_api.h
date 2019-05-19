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

#ifndef MATCH_LAST_API_H
#define MATCH_LAST_API_H

/* The <GtMatchLAST> class, implementing the <GtMatch> interface, is meant to
    store results given in the format as output by LAST. */
typedef struct GtMatchLAST GtMatchLAST;

#include "extended/match_api.h"

/* Creates a new <GtMatch> object meant to store results from the
   LAST software. */
GtMatch*      gt_match_last_new(const char *seqid1,
                                const char *seqid2,
                                GtUword score,
                                GtUword seqno1,
                                GtUword seqno2,
                                GtUword start_seq1,
                                GtUword start_seq2,
                                GtUword end_seq1,
                                GtUword end_seq2,
                                GtMatchDirection dir);

/* Returns the sequence number of the match <ms> in the first sequence set. */
GtUword gt_match_last_get_seqno1(const GtMatchLAST *ml);

/* Returns the sequence number of the match <ms> in the second sequence set. */
GtUword gt_match_last_get_seqno2(const GtMatchLAST *ml);

/* Returns the LAST score of the match <ms>. */
GtUword gt_match_last_get_score(const GtMatchLAST *ml);

#endif
