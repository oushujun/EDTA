/*
  Copyright (c) 2011      Sascha Kastens <mail@skastens.de>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_API_H
#define MATCH_API_H

#include "core/range_api.h"

/* The <GtMatch> interface defines a generic set of functions for a data
   structure designed to hold matches, that is, two similar locations on two
   sequences which are further described by additional data specific to the
   matching engine. Matches have a direction, e.g. direct or reverse
   (palindromic). */
typedef struct GtMatch GtMatch;

typedef enum {
  GT_MATCH_DIRECT,
  GT_MATCH_REVERSE
} GtMatchDirection;

/* Sets the sequence ID of the first sequence involved in the match <match> to
   <seqid>. The string <seqid> must be null-terminated. */
void             gt_match_set_seqid1(GtMatch *match, const char *seqid);
/* Sets the sequence ID of the first sequence involved in the match <match> to
   <seqid>. The string <seqid> needs not be null-terminated, its length is
   given by <len>. */
void             gt_match_set_seqid1_nt(GtMatch *match, const char *seqid,
                                        GtUword len);
/* Sets the sequence ID of the second sequence involved in the match <match> to
   <seqid>. The string <seqid> must be null-terminated. */
void             gt_match_set_seqid2(GtMatch *match, const char *seqid);
/* Sets the sequence ID of the second sequence involved in the match <match> to
   <seqid>. The string <seqid> needs not be null-terminated, its length is
   given by <len>. */
void             gt_match_set_seqid2_nt(GtMatch *match, const char *seqid,
                                        GtUword len);
/* Returns the sequence ID of the first sequence involved in the match
   <match>. */
const char*      gt_match_get_seqid1(const GtMatch *match);
/* Returns the sequence ID of the second sequence involved in the match
   <match>. */
const char*      gt_match_get_seqid2(const GtMatch *match);
/* Sets the range of the first sequence involved in the match <match> to
   <start>-<end>. */
void             gt_match_set_range_seq1(GtMatch *match, GtUword start,
                                         GtUword end);
/* Sets the range of the second sequence involved in the match <match> to
   <start>-<end>. */
void             gt_match_set_range_seq2(GtMatch *match, GtUword start,
                                         GtUword end);
/* Writes the range of the first sequence involved in the match <match> to the
   location pointed to by <range>.
   Note: depending on how the matches were produced the resulting range might
   differ. e.g. Blast hit ranges are 1-Based, not zero based and inclusive i.e.
   <range>.end is the last position that is part of the match. */
void             gt_match_get_range_seq1(const GtMatch *match, GtRange *range);
/* Writes the range of the second sequence involved in the match <match> to the
   location pointed to by <range>.
   Note: depending on how the matches were produced the resulting range might
   differ. e.g. Blast hit ranges are 1-Based, not zero based and inclusive i.e.
   <range>.end is the last position that is part of the match. */
void             gt_match_get_range_seq2(const GtMatch *match, GtRange *range);
/* Returns the match direction of <match>. */
GtMatchDirection gt_match_get_direction(const GtMatch *match);

/* Deletes the match <match>, freeing all associated memory. */
void             gt_match_delete(GtMatch *match);

#endif
