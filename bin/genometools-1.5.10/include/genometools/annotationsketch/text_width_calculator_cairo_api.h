/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef TEXT_WIDTH_CALCULATOR_CAIRO_API_H
#define TEXT_WIDTH_CALCULATOR_CAIRO_API_H

#include <cairo.h>
#include "annotationsketch/style_api.h"
#include "annotationsketch/text_width_calculator_api.h"

/* Implements the GtTextWidthCalculator interface with Cairo as the drawing
   backend. If text width is to be calculated with regard to a specific
   transformation etc. which is in effect in a <cairo_t> and which should be
   used later via a <GtCanvasCairoContext>, create a
   <GtTextWidthCalculatorCairo> object and pass it to the <GtLayout> via
   <gt_layout_new_with_twc()>. */
typedef struct GtTextWidthCalculatorCairo GtTextWidthCalculatorCairo;

/* Creates a new <GtTextWidthCalculatorCairo> object for the given context
   using the text size options given in the <GtStyle>. If the <GtStyle> is NULL,
   the current font settings in the <cairo_t> will be used for all text
   width calculations. */
GtTextWidthCalculator* gt_text_width_calculator_cairo_new(cairo_t*, GtStyle*,
                                                          GtError*);

#endif
