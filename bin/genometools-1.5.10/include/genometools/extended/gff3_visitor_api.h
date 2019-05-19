/*
  Copyright (c) 2006-2008, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008       Center for Bioinformatics, University of Hamburg

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

#ifndef GFF3_VISITOR_API_H
#define GFF3_VISITOR_API_H

#include "core/file_api.h"
#include "extended/node_visitor_api.h"

/* Implements the <GtNodeVisitor> interface with a visitor that produces GFF3
   output. This is a low-level class and it is usually not used directly.
   Normally, a <GtGFF3OutStream> is used to produce GFF3 output. */
typedef struct GtGFF3Visitor GtGFF3Visitor;

/* Create a new <GtNodeVisitor*> which writes the output it produces to the
   given output file pointer <outfp>. If <outfp> is <NULL>, the output is
   written to <stdout>. */
GtNodeVisitor* gt_gff3_visitor_new(GtFile *outfp);
/* Set the width with which the FASTA sequences of <GtSequenceNode>s visited
   by <gff3_visitor> are shown to <fasta_width>.
   Per default, each FASTA entry is shown on a single line. */
void           gt_gff3_visitor_set_fasta_width(GtGFF3Visitor *gff3_visitor,
                                               GtUword fasta_width);
/* Retain the original ID attributes (instead of creating new ones), if
   possible.  Memory consumption for <gff3_visitor> is raised from O(1) to
   O(<input_size>), because bookkeeping of used IDs becomes necessary to avoid
   ID collisions. */
void           gt_gff3_visitor_retain_id_attributes(GtGFF3Visitor
                                                    *gff3_visitor);

#endif
