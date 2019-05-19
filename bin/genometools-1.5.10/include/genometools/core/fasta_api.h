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

#ifndef FASTA_API_H
#define FASTA_API_H

#include "core/file_api.h"
#include "core/str_api.h"
#include "core/types_api.h"

# define GT_FASTA_DEFAULT_WIDTH ((GtUword) 80)

/* FASTA module */

/* Print a fasta entry with optional <description> and mandatory <sequence> to
   <outfp>. If <width> is != 0 the sequence is formatted accordingly. */
void gt_fasta_show_entry(const char *description, const char *sequence,
                         GtUword sequence_length, GtUword width, GtFile *outfp);

/* Print a fasta entry with optional <description> and mandatory <sequence> to
   <outfp>. If <width> is != 0 the sequence is formatted accordingly.
   Will print at most <sequence_length> characters from <sequence> if no
   '\0'-byte is encountered. <description> and <description_length> are handled
   accordingly if present. */
void gt_fasta_show_entry_nt(const char *description, GtUword description_length,
                            const char *sequence, GtUword sequence_length,
                            GtUword width, GtFile *outfp);

/* Append a fasta entry with optional <description> and mandatory <sequence> to
   <outstr>. If <width> is != 0 the sequence is formatted accordingly. */
void gt_fasta_show_entry_str(const char *description, const char *sequence,
                             GtUword sequence_length, GtUword width,
                             GtStr *outstr);

/* Print a fasta entry with optional <description> and mandatory <sequence> to
   <outstr>. If <width> is != 0 the sequence is formatted accordingly.
   Will print at most <sequence_length> characters from <sequence> if no
   '\0'-byte is encountered. <description> and <description_length> are handled
   accordingly if present. */
void gt_fasta_show_entry_nt_str(const char *description,
                                GtUword description_length,
                                const char *sequence,
                                GtUword sequence_length,
                                GtUword width,
                                GtStr *outstr);
#endif
