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

#ifndef CODON_ITERATOR_SIMPLE_API_H
#define CODON_ITERATOR_SIMPLE_API_H

#include "core/codon_iterator_api.h"
#include "core/error_api.h"

typedef struct GtCodonIteratorSimple GtCodonIteratorSimple;

/* Creates a new <GtCodonIterator> traversing <seq> over a length of <len>.
   If an error occurs, NULL is returned and <err> is set accordingly. */
GtCodonIterator*            gt_codon_iterator_simple_new(const char *seq,
                                                         GtUword len,
                                                         GtError *err);

const GtCodonIteratorClass* gt_codon_iterator_simple_class(void);
int                         gt_codon_iterator_simple_unit_test(GtError *err);

#endif
