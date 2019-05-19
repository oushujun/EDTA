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

#ifndef CANVAS_CAIRO_CONTEXT_API_H
#define CANVAS_CAIRO_CONTEXT_API_H

#include <cairo.h>
#include "annotationsketch/canvas_api.h"
#include "annotationsketch/image_info_api.h"
#include "annotationsketch/style_api.h"

/* Implements the <GtCanvas> interface using a Cairo context (<cairo_t>)
   as input. This Canvas uses the <GtGraphicsCairo> class.

   Drawing to a <cairo_t> allows the use of the  __AnnotationSketch__ engine
   in any Cairo-based graphical application. */
typedef struct GtCanvasCairoContext GtCanvasCairoContext;

/* Create a new <GtCanvas> object tied to the cairo_t <context>, <width> and
   <height> using the given <style>. The optional <image_info> is
   filled when the created Canvas object is used to render a <GtDiagram> object.
   <offsetpos> determines where to start drawing on the surface. */
GtCanvas* gt_canvas_cairo_context_new(GtStyle *style, cairo_t *context,
                                      double offsetpos,
                                      GtUword width,
                                      GtUword height,
                                      GtImageInfo *image_info,
                                      GtError *err);
#endif
