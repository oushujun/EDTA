/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
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

#ifndef GFF3_IN_STREAM_API_H
#define GFF3_IN_STREAM_API_H

#include <stdio.h>
#include "core/str_array_api.h"
#include "extended/node_stream_api.h"
#include "extended/type_checker_api.h"

/* Implements the <GtNodeStream> interface. A <GtGFF3InStream> parses GFF3 files
   and returns them as a stream of <GtGenomeNode> objects. */
typedef struct GtGFF3InStream GtGFF3InStream;

/* Return a <GtGFF3InStream> object which subsequently reads the <num_of_files>
   many GFF3 files denoted in <filenames>. The GFF3 files do not have to be
   sorted. If <num_of_files> is 0 or a file name is "-", it is read from
   <stdin>. The memory footprint is O(file size) in the worst-case. */
GtNodeStream* gt_gff3_in_stream_new_unsorted(int num_of_files,
                                             const char **filenames);
/* Create a <GtGFF3InStream*> which reads the sorted GFF3 file denoted by
   <filename>. If filename is <NULL>, it is read from <stdin>.
   The memory footprint is O(1) on average. */
GtNodeStream* gt_gff3_in_stream_new_sorted(const char *filename);
/* Make sure all ID attributes which are parsed by <gff3_in_stream> are correct.
   Increases the memory footprint to O(file size). */
void          gt_gff3_in_stream_check_id_attributes(GtGFF3InStream
                                                               *gff3_in_stream);
/* Enable tidy mode for <gff3_in_stream>. That is, the GFF3 parser tries to tidy
   up features which would normally lead to an error. */
void          gt_gff3_in_stream_enable_tidy_mode(GtGFF3InStream
                                                               *gff3_in_stream);
/* Enable strict mode for <gff3_in_stream>. */
void          gt_gff3_in_stream_enable_strict_mode(GtGFF3InStream
                                                               *gff3_in_stream);
/* Show progress bar on <stdout> to convey the progress of parsing the GFF3
   files underlying <gff3_in_stream>. */
void          gt_gff3_in_stream_show_progress_bar(GtGFF3InStream
                                                               *gff3_in_stream);
/* Returns a <GtStrArray*> which contains all type names in alphabetical order
   which have been parsed by <gff3_in_stream>. The caller is responsible to
   free it! */
GtStrArray*   gt_gff3_in_stream_get_used_types(GtNodeStream *gff3_in_stream);
/* Sets <type_checker> to be the type checker used in <gff3_in_stream>. That
   is, it will be queried when the validity of SO types is to be determined. */
void          gt_gff3_in_stream_set_type_checker(GtNodeStream *gff3_in_stream,
                                                 GtTypeChecker *type_checker);

#endif
