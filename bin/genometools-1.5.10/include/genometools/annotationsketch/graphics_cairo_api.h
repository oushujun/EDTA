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

#ifndef GRAPHICS_CAIRO_API_H
#define GRAPHICS_CAIRO_API_H

#include <cairo.h>
#include "annotationsketch/graphics_api.h"

/* Implements the <GtGraphics> interface.
   This implementation uses the Cairo 2D vector graphics library as a
   drawing back-end. */
typedef struct GtGraphicsCairo GtGraphicsCairo;

/* Creates a new <GtGraphics> object using the Cairo backend. The object
   is meant for writing a new image of width <width> and height <height>
   to a file or stream. Use <type> to define the output format.  */
GtGraphics*  gt_graphics_cairo_new(GtGraphicsOutType type, unsigned int width,
                                   unsigned int height);
/* Creates a new <GtGraphics> object using the Cairo backend. The object
   is meant for writing on an existing cairo_t <context> within the boundaries
   of width <width> and height <height>.  */
GtGraphics*  gt_graphics_cairo_new_from_context(cairo_t *context,
                                                unsigned int width,
                                                unsigned int height);
/* Draws a curve in <gg> at the position <x>,<y> for <ndata> data points as
   given in <data>. The data points must be in the range <valrange> and the
   resulting graph has the height <height> in type-dependent units
   (e.g. pixels). */
void         gt_graphics_cairo_draw_curve_data(GtGraphics *gg,
                                               double x, double y,
                                               GtColor color, double data[],
                                               GtUword ndata, GtRange valrange,
                                               GtUword height);
#endif
