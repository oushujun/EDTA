/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef TRANSLATOR_API_H
#define TRANSLATOR_API_H

#include "core/codon_iterator_api.h"
#include "core/error_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "core/trans_table_api.h"

typedef enum {
  GT_TRANSLATOR_OK    =  0,
  GT_TRANSLATOR_END   = -1,
  GT_TRANSLATOR_ERROR = -2
} GtTranslatorStatus;

/* The <GtTranslator> can be used to  produce 3-frame translations of DNA
   sequences via an iterator interface. */
typedef struct GtTranslator GtTranslator;

/* Creates a new <GtTranslator>, starting its translation at the current
   position of <ci>. The current reading frame is also taken from the state of
   <ci>. The translation table <tt> is used. */
GtTranslator*      gt_translator_new_with_table(GtTransTable *tt,
                                                GtCodonIterator *ci);

/* Creates a new <GtTranslator>, starting its translation at the current
   position of <ci>. The current reading frame is also taken from the state of
   <ci>. The standard translation table is used. */
GtTranslator*      gt_translator_new(GtCodonIterator *ci);

/* Reinitializes <translator> with the position and frame status as given in
   <ci>. */
void               gt_translator_set_codon_iterator(GtTranslator *translator,
                                                    GtCodonIterator *ci);

/* Selects the translation scheme in <translator> to the one identified by
   translation table <tt>. */
void               gt_translator_set_translation_table(GtTranslator *translator,
                                                       GtTransTable *tt);

/* Returns the translation of the next codon. The currently translated
   character is put in <translated> while the current reading frame is put in
   <frame>.
   Returns GT_TRANSLATOR_ERROR if an error occurred, see <err> for details.
   If the end of the sequence region to translate has been reached,
   GT_TRANSLATOR_END is returned.
   Otherwise, GT_TRANSLATOR_OK (equal to 0) is returned. */
GtTranslatorStatus gt_translator_next(GtTranslator *translator,
                                      char *translated,
                                      unsigned int *frame,
                                      GtError *err);

/* Moves the <translator> to the beginning of the first codon in <dnaseq> (of
   length <dnalen>) which is a start codon according to the selected translation
   scheme in <translator>.
   The offset is written to the location pointed to by <pos>.
   Returns GT_TRANSLATOR_ERROR if an error occurred, see <err> for details.
   If the end of the sequence region to scan has been reached without finding a
   start codon, GT_TRANSLATOR_END is returned.
   Otherwise, GT_TRANSLATOR_OK (equal to 0) is returned. */
GtTranslatorStatus gt_translator_find_startcodon(GtTranslator *translator,
                                                 GtUword *pos,
                                                 GtError *err);

/* Moves the <translator> to the beginning of the first codon in <dnaseq> (of
   length <dnalen>) which is a stop codon according to the selected translation
   scheme in <translator>.
   The offset is written to the location pointed to by <pos>.
   Returns GT_TRANSLATOR_ERROR if an error occurred, see <err> for details.
   If the end of the sequence region to scan has been reached without finding a
   stop codon, GT_TRANSLATOR_END is returned.
   Otherwise, GT_TRANSLATOR_OK (equal to 0) is returned. */
GtTranslatorStatus gt_translator_find_stopcodon(GtTranslator *translator,
                                                GtUword *pos,
                                                GtError *err);

/* Moves the <translator> to the beginning of the first codon in <dnaseq> (of
   length <dnalen>) which belongs to the set of codons specified in <codons>.
   The offset is written to the location pointed to by <pos>.
   Returns GT_TRANSLATOR_ERROR if an error occurred, see <err> for details.
   If the end of the sequence region to scan has been reached without finding
   one of the codons, GT_TRANSLATOR_END is returned.
   Otherwise, GT_TRANSLATOR_OK (equal to 0) is returned. */
GtTranslatorStatus gt_translator_find_codon(GtTranslator *translator,
                                            GtStrArray *codons,
                                            GtUword *pos,
                                            GtError *err);

/* Delete <translator>. */
void               gt_translator_delete(GtTranslator *translator);

#endif
