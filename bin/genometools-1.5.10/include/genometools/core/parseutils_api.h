/*
  Copyright (c) 2006-2012 Gordon Gremme <gordon@gremme.org>
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

#ifndef PARSEUTILS_API_H
#define PARSEUTILS_API_H

#include "core/deprecated_api.h"
#include "core/error_api.h"
#include "core/phase_api.h"
#include "core/range_api.h"
#include "core/strand_api.h"

/* Parseutils module */

/* Parse integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_int(int *out, const char *nptr);

/* Parse unsigned integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_uint(unsigned int *out, const char *nptr);

/* Parse long from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
GT_DEPRECATED("use gt_parse_word() instead")
int gt_parse_long(GtWord *out, const char *nptr);
/* Parse GtWord from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_word(GtWord *out, const char *nptr);

/* Parse ulong from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
GT_DEPRECATED("use gt_parse_uword() instead")
int gt_parse_ulong(GtUword *out, const char *nptr);
/* Parse GtUword from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_uword(GtUword *out, const char *nptr);

/* Parse double from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_double(double *out, const char *nptr);

/* Parse a range given by <start> and <end>, writing the result into <rng>.
   Enforces that <start> is smaller or equal than <end>. Give <filename> and
   <line_number> for error reporting. Returns 0 upon success and -1 upon
   failure. */
int gt_parse_range(GtRange *rng, const char *start, const char *end,
                   unsigned int line_number, const char *filename, GtError*);

/* Like <gt_parse_range>, but issues a warning if <start> is larger then <end>
   and swaps both values. It also issues a warning, if <start> and/or <end> is
   not-positive and sets the corresponding value to 1. */
int gt_parse_range_tidy(GtRange *rng, const char *start, const char *end,
                        unsigned int line_number, const char *filename,
                        GtError*);

#endif
