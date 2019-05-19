/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef DISC_DISTRI_API_H
#define DISC_DISTRI_API_H

#include "core/error_api.h"
#include "core/file_api.h"
#include "core/types_api.h"

/* The <GtDiscDistri> class represents a discrete distribution of integer
   values. */
typedef struct GtDiscDistri GtDiscDistri;

/* Callback function called during iteration over each item of the
   distribution, where <key> is the counted value and <value> is the count. */
typedef void (*GtDiscDistriIterFunc)(GtUword key,
                                     GtUint64 value,
                                     void *data);

/* Creates a new, empty <GtDiscDistri>. */
GtDiscDistri* gt_disc_distri_new(void);
/* Adds one count of <key> to <d>. */
void          gt_disc_distri_add(GtDiscDistri *d, GtUword key);
/* Adds <occurrences> counts of <key> to <d>. */
void          gt_disc_distri_add_multi(GtDiscDistri *d, GtUword key,
                                       GtUint64 occurrences);
/* Return the current count of <key> as stored in <d>. */
GtUint64      gt_disc_distri_get(const GtDiscDistri *d, GtUword key);
/* Prints the current state of <d> to <outfile>. If <outfp> is NULL,
   stdout will be used for output. */
void          gt_disc_distri_show(const GtDiscDistri *d, GtFile *outfp);
/* Iterate over all non-empty entries in <d>, calling <func> for each one,
   from the smallest to the largest key. The <data> pointer can be used to pass
   arbitrary data to <func>. */
void          gt_disc_distri_foreach(const GtDiscDistri *d,
                                     GtDiscDistriIterFunc func,
                                     void *data);
/* Same as foreach, but from the longest to the smallest key. */
void          gt_disc_distri_foreach_in_reverse_order(const GtDiscDistri *d,
                                                      GtDiscDistriIterFunc func,
                                                      void *data);
int           gt_disc_distri_unit_test(GtError*);
void          gt_disc_distri_delete(GtDiscDistri*);

#endif
