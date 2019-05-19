/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef CODON_ITERATOR_API_H
#define CODON_ITERATOR_API_H

#include "core/error_api.h"
#include "core/types_api.h"

typedef enum {
  GT_CODON_ITERATOR_OK    =  0,
  GT_CODON_ITERATOR_END   = -1,
  GT_CODON_ITERATOR_ERROR = -2
} GtCodonIteratorStatus;

/* The <GtCodonIterator> interface. */
typedef struct GtCodonIterator GtCodonIterator;
typedef struct GtCodonIteratorClass GtCodonIteratorClass;

/* Return the current reading offset of <codin_iterator>, starting from the
   position in the sequence given at iterator instantiation time. */
GtUword          gt_codon_iterator_current_position(GtCodonIterator
                                                          *codon_iterator);
/* Return the length of the substring to scan, given at instantiation time. */
GtUword          gt_codon_iterator_length(GtCodonIterator
                                                *codon_iterator);
/* Rewind the <codon_iterator> to point again to the position in the sequence
   given at iterator instantiation time. */
void                   gt_codon_iterator_rewind(GtCodonIterator
                                                *codon_iterator);
/* Sets the values of <n1>, <n2> and <n3> to the codon beginning at the current
   reading position of <codon_iterator> and then advances the reading position
   by one. The current reading frame shift (0, 1 or 2) is for the current codon
   is written to the position pointed to by <frame>.
   This function returns one of three status codes:
   GT_CODON_ITERATOR_OK    : a codon was read successfully,
   GT_CODON_ITERATOR_END   : no codon was read because the end of the scan
                             region has been reached,
   GT_CODON_ITERATOR_ERROR : no codon was read because an error occurred during
                             sequence access. See <err> for details. */
GtCodonIteratorStatus  gt_codon_iterator_next(GtCodonIterator *codon_iterator,
                                              char *n1, char *n2, char *n3,
                                              unsigned int *frame,
                                              GtError *err);
/* Delete <codon_iterator>. */
void                   gt_codon_iterator_delete(GtCodonIterator
                                                *codon_iterator);

#endif
