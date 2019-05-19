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

#ifndef TRANS_TABLE_API_H
#define TRANS_TABLE_API_H

#include "core/error_api.h"
#include "core/str_array_api.h"

typedef struct GtTransTable GtTransTable;

/* The number of the standard translation scheme. */
#define GT_STANDARD_TRANSLATION_SCHEME 1U

/* Returns a <GtStrArray> of translation scheme descriptions, each of the
   format "%d: %s" where the number is the translation scheme number (usable in
   <gt_translator_set_translation_scheme()> and the string is the scheme
   name. */
GtStrArray*   gt_trans_table_get_scheme_descriptions(void);

/* Returns a translation table as given by <scheme> which refers to the numbers
   as reported by <gt_translator_get_translation_table_descriptions()> or the
   list given at the NCBI web site
   __http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi__.
   Returns NULL if an error occurred, see <err> for details. */
GtTransTable* gt_trans_table_new(unsigned int scheme, GtError *err);

/* Returns the standard translation table. */
GtTransTable* gt_trans_table_new_standard(GtError *err);

/* Returns the description of <tt>. */
const char*   gt_trans_table_description(const GtTransTable *tt);

/* Writes the translation for the codon <c1>,<c2>,<c3> to the position pointed
   to by <amino>. The current translation scheme set in <translator> is used.
   Returns a negative value if an error occurred, see <err> for details.
   Otherwise, 0 is returned. */
int           gt_trans_table_translate_codon(const GtTransTable *tt,
                                             char c1, char c2, char c3,
                                             char *amino, GtError *err);

bool          gt_trans_table_is_start_codon(const GtTransTable *tt,
                                            char c1, char c2, char c3);

bool          gt_trans_table_is_stop_codon(const GtTransTable *tt,
                                           char c1, char c2, char c3);

/* Deletes <tt>. */
void          gt_trans_table_delete(GtTransTable *tt);

#endif
