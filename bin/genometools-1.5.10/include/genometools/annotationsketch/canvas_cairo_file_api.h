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

#ifndef CANVAS_CAIRO_FILE_API_H
#define CANVAS_CAIRO_FILE_API_H

#include "annotationsketch/canvas_api.h"
#include "annotationsketch/graphics_api.h"
#include "annotationsketch/image_info_api.h"
#include "annotationsketch/style_api.h"

/* Implements the <GtCanvas> interface.
   This Canvas uses the <GtGraphicsCairo> class.  */
typedef struct GtCanvasCairoFile GtCanvasCairoFile;

/* Create a new <GtCanvasCairoFile> object with given <output_type> and
   <width> using the configuration given in <style>. The optional <image_info>
   is filled when the created object is used to render a <GtDiagram> object.
   Possible <GtGraphicsOutType> values are <GRAPHICS_PNG>, <GRAPHICS_PS>,
   <GRAPHICS_PDF> and <GRAPHICS_SVG>. Dependent on the local Cairo installation,
   not all of them may be available. */
GtCanvas* gt_canvas_cairo_file_new(GtStyle *style,
                                   GtGraphicsOutType output_type,
                                   GtUword width,
                                   GtUword height,
                                   GtImageInfo *image_info,
                                   GtError *err);
/* Write rendered <canvas> to the file with name <filename>. If this
   method returns a value other than 0, check <err> for an error message. */
int     gt_canvas_cairo_file_to_file(GtCanvasCairoFile *canvas,
                                     const char *filename, GtError *err);
/* Append rendered <canvas> image data to given <stream>. */
int     gt_canvas_cairo_file_to_stream(GtCanvasCairoFile *canvas,
                                       GtStr *stream);

#endif
