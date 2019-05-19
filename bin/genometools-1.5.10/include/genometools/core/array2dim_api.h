/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef ARRAY2DIM_API_H
#define ARRAY2DIM_API_H

#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/types_api.h"

/* Array2dim module */

/* Allocates a new 2-dimensional array with dimensions <ROWS> x <COLUMNS> and
   assigns a pointer to the newly allocated space to <ARRAY2DIM>.
   The size of each element is determined automatically from the type of the
   <ARRAY2DIM> pointer. */
#define gt_array2dim_malloc(ARRAY2DIM, ROWS, COLUMNS)                          \
        {                                                                      \
          GtUword gt_a2d_i;                                                    \
          ARRAY2DIM = gt_malloc(sizeof *ARRAY2DIM * (ROWS));                   \
          (ARRAY2DIM)[0] = gt_malloc(sizeof **ARRAY2DIM * (ROWS) * (COLUMNS)); \
          for (gt_a2d_i = 1UL; gt_a2d_i < (GtUword) (ROWS); gt_a2d_i++)        \
            (ARRAY2DIM)[gt_a2d_i] = (ARRAY2DIM)[gt_a2d_i-1] + (COLUMNS);       \
        }

/* Allocates a new 2-dimensional array with dimensions <ROWS> x <COLUMNS> and
   assigns a pointer to the newly allocated space to <ARRAY2DIM>.
   The allocated space is initialized to be filled with zeroes.
   The size of each element is determined automatically from the type of the
   <ARRAY2DIM> pointer. */
#define gt_array2dim_calloc(ARRAY2DIM, ROWS, COLUMNS)                         \
        {                                                                     \
          GtUword gt_a2d_i;                                                   \
          ARRAY2DIM = gt_malloc(sizeof *ARRAY2DIM * (ROWS));                  \
          (ARRAY2DIM)[0] = gt_calloc((size_t) ((ROWS) * (COLUMNS)),           \
            sizeof **ARRAY2DIM);                                              \
          for (gt_a2d_i = 1UL; gt_a2d_i < (GtUword) (ROWS); gt_a2d_i++)       \
            (ARRAY2DIM)[gt_a2d_i] = (ARRAY2DIM)[gt_a2d_i-1] + (COLUMNS);      \
        }

/* An example for usage of the <Array2dim> module. */
int     gt_array2dim_example(GtError*);

/* Frees the space allocated for the 2-dimensional array pointed to by
   <ARRAY2DIM>. */
#define gt_array2dim_delete(ARRAY2DIM) \
        gt_free((ARRAY2DIM)[0]);       \
        gt_free(ARRAY2DIM);

#endif
