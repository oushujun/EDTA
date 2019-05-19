/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TEXT_WIDTH_CALCULATOR_API_H
#define TEXT_WIDTH_CALCULATOR_API_H

#include "core/error_api.h"

/* The GtTextWidthCalculator interface answers queries w.r.t.
   text width in a specific drawing backend. This interface is needed to do
   proper line breaking in a <GtLayout> even if there is no <GtCanvas> or
   <GtGraphics> created yet. */
typedef struct GtTextWidthCalculator GtTextWidthCalculator;

/* Increases the reference count of the <GtTextWidthCalculator>. */
GtTextWidthCalculator* gt_text_width_calculator_ref(GtTextWidthCalculator*);
/* Requests the width of <text> from the <GtTextWidthCalculator>.
   If the returned value is negative, an error occurred. Otherwise,
   a positive double value is returned. */
double                 gt_text_width_calculator_get_text_width(
                                                    GtTextWidthCalculator*,
                                                    const char *text,
                                                    GtError *err);
/* Deletes a <GtTextWidthCalculator> instance. */
void                   gt_text_width_calculator_delete(GtTextWidthCalculator*);

#endif
