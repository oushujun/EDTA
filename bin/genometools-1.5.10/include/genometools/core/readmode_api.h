/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef READMODE_API_H
#define READMODE_API_H

#include "core/error_api.h"

enum GtReadmode {
  GT_READMODE_FORWARD = 0,
  GT_READMODE_REVERSE,
  GT_READMODE_COMPL,
  GT_READMODE_REVCOMPL
};

/* This enum type defines the possible reamodes, namely <GT_READMODE_FORWARD>,
   <GT_READMODE_REVERSE>, <GT_READMODE_COMPL>, and <GT_READMODE_REVCOMPL>. */
typedef enum GtReadmode GtReadmode;

#define gt_readmode_invert(RM) \
        RM = ((GtReadmode) (3 - (int) (RM)))

/* Returns the descriptive string for <readmode>. */
const char* gt_readmode_show(GtReadmode readmode);
/* Returns the <GtReadmode> for the description <string>, which must be one
   of "fwd","rev","cpl" or "rcl". If <string> does not equal any of them,
   -1 is returned and <err> is set accordingly. */
int         gt_readmode_parse(const char *string, GtError *err);

/* invert the direction of a readmode: fwd => rev, rev => fwd,
                                       cpl => rcl, rcl => cpl */
GtReadmode gt_readmode_inverse_dir(GtReadmode readmode);

#endif
