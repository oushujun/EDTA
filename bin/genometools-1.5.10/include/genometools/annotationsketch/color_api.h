/*
  Copyright (c) 2007 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef COLOR_API_H
#define COLOR_API_H

#include <stdbool.h>

/* The <GtColor> class holds a RGB color definition. */
typedef struct GtColor GtColor;

struct GtColor {
  double red, green, blue, alpha;
};

/* Create a new <GtColor> object with the color given by the <red>, <green>,
   and <blue> arguments. The value for each color channel must be between 0
   and 1. */
GtColor* gt_color_new(double red, double green, double blue, double alpha);
/* Change the color of the <color> object to the color given by the <red>,
   <green>, and <blue> arguments. The value for each color channel must be
   between 0 and 1.*/
void      gt_color_set(GtColor *color, double red, double green, double blue,
                       double alpha);
/* Returns <true> if the colors <c1> and <c2> are equal. */
bool      gt_color_equals(const GtColor *c1, const GtColor *c2);
/* Delete the <color> object. */
void      gt_color_delete(GtColor *color);

#endif
