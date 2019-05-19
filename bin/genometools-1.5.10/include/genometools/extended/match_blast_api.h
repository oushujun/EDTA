/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
  Copyright (c)      2016 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2011-2016 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_BLAST_API_H
#define MATCH_BLAST_API_H

typedef struct GtMatchBlast GtMatchBlast;

#include "extended/match_api.h"

/* Creates a new <GtMatch> object meant to store results in the BLAST
   format. That is, it stores double values <evalue> for match E-values,
   <bitscore>s and the alignment length <ali_l> in addition to the generic
   match contents <seqid1>, <seqid2>, <start_seq1>, <start_seq2>, <end_seq1>,
   and <end_seq2>. */
GtMatch* gt_match_blast_new(char *seqid1,
                            char *seqid2,
                            GtUword start_seq1,
                            GtUword start_seq2,
                            GtUword end_seq1,
                            GtUword end_seq2,
                            double evalue,
                            float bitscore,
                            GtUword ali_l,
                            double similarity,
                            GtMatchDirection dir);

/* Creates a new <GtMatch> object meant to store results in the BLAST
   format. That is, it stores double values <evalue> for match E-values,
   <bitscore>s and the alignment length <ali_l> in addition to the generic
   match contents <seqid1>, <seqid2>, <start_seq1>, <start_seq2>, <end_seq1>,
   and <end_seq2>. In addition to <gt_match_blast_new> it also stores
   the number of mismatches and the number of gap */

GtMatch* gt_match_blast_new_extended(char *seqid1,
                                     char *seqid2,
                                     GtUword start_seq1,
                                     GtUword end_seq1,
                                     GtUword start_seq2,
                                     GtUword end_seq2,
                                     double evalue,
                                     float bitscore,
                                     GtUword length,
                                     double similarity,
                                     GtUword mm_num,
                                     GtUword gap_open_num,
                                     GtMatchDirection dir);

/* Sets <evalue> to be the E-value in <mb>. */
void gt_match_blast_set_evalue(GtMatchBlast *mb, double evalue);

/* Sets <bits> to be the bit-score in <mb>. */
void gt_match_blast_set_bitscore(GtMatchBlast *mb, float bits);

/* Sets <length> to be the alignment length in <mb>. */
void gt_match_blast_set_align_length(GtMatchBlast *mb, GtUword length);

/* Sets <similarity> to be the match similarity in <mb>. */
void gt_match_blast_set_similarity(GtMatchBlast *mb, double similarity);

/* Sets <num to be the number of mismatches in <mb>. */
void gt_match_blast_set_mismatches(GtMatchBlast *mb, GtUword mm_num);

/* Sets <num> to be the number of gap openings in <mb>. */
void gt_match_blast_set_gapopen(GtMatchBlast *mb, GtUword gap_open_num);

/* Returns the E-value stored in <mb>. */
double gt_match_blast_get_evalue(GtMatchBlast *mb);

/* Returns the bit-score value stored in <mb>. */
float gt_match_blast_get_bitscore(GtMatchBlast *mb);

/* Returns the alignment length stored in <mb>. */
GtUword gt_match_blast_get_align_length(GtMatchBlast *mb);

/* Returns the alignment similarity stored in <mb>. */
double gt_match_blast_get_similarity(GtMatchBlast *mb);

/* Returns the number of mismatches stored in <mb>. */
GtUword gt_match_blast_get_mismatches(GtMatchBlast *mb);

/* Returns the number of gap openings stored in <mb>. */
GtUword gt_match_blast_get_gapopen(GtMatchBlast *mb);
#endif
