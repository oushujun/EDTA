/*
  Copyright (c) 2009-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQ_ITERATOR_FASTQ_API_H
#define SEQ_ITERATOR_FASTQ_API_H

#include "core/error_api.h"
#include "core/seq_iterator_api.h"
#include "core/str_array_api.h"

typedef struct GtSeqIteratorFastQ GtSeqIteratorFastQ;

/* Create a new <GtSeqIteratorFastQ> for all sequence files in <filenametab>. */
GtSeqIterator* gt_seq_iterator_fastq_new(const GtStrArray *filenametab,
                                         GtError *err);

/* Create a new <GtSeqIteratorFastQ> for all sequence files in <filenametab>
   containing color space reads. */
GtSeqIterator* gt_seq_iterator_fastq_new_colorspace(const GtStrArray
                                                    *filenametab,
                                                    GtError *err);

/* Returns the number of the file in the file name array which <seqit> is
   currently reading.  */
GtUword  gt_seq_iterator_fastq_get_file_index(GtSeqIteratorFastQ *seqit);

/* Disable checking if quality description is equal to read description in
   <seqit> (it should be, but it is not in output of some tools, e.g. Coral). */
void           gt_seq_iterator_fastq_relax_check_of_quality_description(
                                                    GtSeqIteratorFastQ *seqit);

const GtSeqIteratorClass* gt_seq_iterator_fastq_class(void);

#endif
