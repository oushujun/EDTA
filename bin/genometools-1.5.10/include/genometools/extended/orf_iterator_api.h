/*
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef ORF_ITERATOR_API_H
#define ORF_ITERATOR_API_H

#include "core/codon_iterator_api.h"
#include "core/translator_api.h"
#include "core/range_api.h"
#include "core/error_api.h"

typedef enum {
  GT_ORF_ITERATOR_OK    =  0,
  GT_ORF_ITERATOR_END   = -1,
  GT_ORF_ITERATOR_ERROR = -2
} GtORFIteratorStatus;

/* The <GtORFIterator> class is used to enumerate open reading frames (ORFs)
   in codon streams as delivered by a <GtCodonIterator> using a translation
   process as defined by a <GtTranslator>. */
typedef struct GtORFIterator GtORFIterator;

/* Return a new <GtORFIterator*> using the codons delivered by <ci> and
   translated by <translator>. */
GtORFIterator*      gt_orf_iterator_new(GtCodonIterator *ci,
                                        GtTranslator *translator);

/* Sets the values of <orf_rng.start>, <orf_rng.end> and <orf_frame> to the
   current reading position of <ci> if a START/STOP codon is found. The frame
   in which the ORF is located is written to the position pointed to by
   <orf_frame>. This function returns one of three status codes:
   <GT_ORF_ITERATOR_OK> if an ORF was detected successfully(START/STOP AA pair),
   <GT_ORF_ITERATOR_END> if no ORF was detected because the end of the scan
                           region has been reached, or
   <GT_ORF_ITERATOR_ERROR> if no ORF was detected because an error occurred
                           during sequence access. See <err> for details. */
GtORFIteratorStatus gt_orf_iterator_next(GtORFIterator *orf_iterator,
                                         GtRange *orf_rng,
                                         unsigned int *orf_frame,
                                         GtError *err);

/* Delete <orf_iterator> and frees all associated memory. */
void                gt_orf_iterator_delete(GtORFIterator *orf_iterator);

#endif
