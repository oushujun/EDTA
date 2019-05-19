/*
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef BITTAB_API_H
#define BITTAB_API_H

#include <stdbool.h>
#include <stdio.h>
#include "core/array_api.h"

/* Implements arbitrary-length bit arrays and various operations on them. */
typedef struct GtBittab GtBittab;

/* Return a new <GtBittab> of length <num_of_bits>, initialised to 0.
   <num_of_bits> has to be > 0 */
GtBittab* gt_bittab_new(GtUword num_of_bits);

/* Set bit <i> in <bittab> to 1. */
void      gt_bittab_set_bit(GtBittab *bittab, GtUword i);

/* Set bit <i> in <bittab> to 0. */
void      gt_bittab_unset_bit(GtBittab *bittab, GtUword i);

/* Set <bittab_a> to be the complement of <bittab_b>. */
void      gt_bittab_complement(GtBittab *bittab_a, const GtBittab *bittab_b);

/* Set <bittab_a> to be equal to <bittab_b>. */
void      gt_bittab_equal(GtBittab *bittab_a, const GtBittab *bittab_b);

/* Set <bittab_a> to be the bitwise AND of <bittab_b> and <bittab_c>. */
void      gt_bittab_and(GtBittab *bittab_a, const GtBittab *bittab_b,
                        const GtBittab *bittab_c);

/* Set <bittab_a> to be the bitwise OR of <bittab_b> and <bittab_c>. */
void      gt_bittab_or(GtBittab *bittab_a, const GtBittab *bittab_b,
                       const GtBittab *bittab_c);

/* Set <bittab_a> to be <bittab_b> NAND <bittab_c>. */
void      gt_bittab_nand(GtBittab *bittab_a, const GtBittab *bittab_b,
                         const GtBittab *bittab_c);

/* Set <bittab_a> to be the bitwise AND of <bittab_a> and <bittab_b>. */
void      gt_bittab_and_equal(GtBittab *bittab_a, const GtBittab *bittab_b);

/* Set <bittab_a> to be the bitwise OR of <bittab_a> and <bittab_b>. */
void      gt_bittab_or_equal(GtBittab *bittab_a, const GtBittab *bittab_b);

/* Shift <bittab> by one position to the left. */
void      gt_bittab_shift_left_equal(GtBittab *bittab);

/* Shift <bittab> by one position to the right. */
void      gt_bittab_shift_right_equal(GtBittab *bittab);

/* Set all bits in <bittab> to 0. */
void      gt_bittab_unset(GtBittab *bittab);

/* Output a representation of <bittab> to <fp>. */
void      gt_bittab_show(const GtBittab *bittab, FILE *fp);

/* Fill <array> with the indices of all set bits in <bittab>. */
void      gt_bittab_get_all_bitnums(const GtBittab *bittab, GtArray *array);

/* Return <true> if bit <i> is set in <bittab>. */
bool      gt_bittab_bit_is_set(const GtBittab *bittab, GtUword i);

/* Return <true> if <bittab_a> and <bittab_b> are identical. */
bool      gt_bittab_cmp(const GtBittab *bittab_a, const GtBittab *bittab_b);

/* Return the index of the first set bit in <bittab>. */
GtUword   gt_bittab_get_first_bitnum(const GtBittab *bittab);

/* Return the index of the last set bit in <bittab>. */
GtUword   gt_bittab_get_last_bitnum(const GtBittab *bittab);

/* Return the index of the next set bit in <bittab> with an index greater
   than <i>. */
GtUword   gt_bittab_get_next_bitnum(const GtBittab *bittab, GtUword i);

/* Return the number of set bits in <bittab>. */
GtUword   gt_bittab_count_set_bits(const GtBittab *bittab);

/* Return the total number of bits of <bittab>. */
GtUword   gt_bittab_size(GtBittab *bittab);

/* Delete <bittab>. */
void      gt_bittab_delete(GtBittab *bittab);

#endif
