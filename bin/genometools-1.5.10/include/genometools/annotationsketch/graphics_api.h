/*
  Copyright (c) 2007-2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>,
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GRAPHICS_API_H
#define GRAPHICS_API_H

#include "annotationsketch/color_api.h"
#include "core/error_api.h"
#include "core/range_api.h"
#include "core/str_api.h"

typedef enum {
  GT_GRAPHICS_PDF,
  GT_GRAPHICS_PNG,
  GT_GRAPHICS_PS,
  GT_GRAPHICS_SVG
} GtGraphicsOutType;

typedef enum
{
  ARROW_LEFT,
  ARROW_RIGHT,
  ARROW_BOTH,
  ARROW_NONE
} ArrowStatus;

typedef enum
{
  SLANT_NORMAL,
  SLANT_ITALIC
} FontSlant;

typedef enum
{
  WEIGHT_NORMAL,
  WEIGHT_BOLD
} FontWeight;

/* The <GtGraphics> interface acts as a low-level abstraction of a drawing
   surface. It is used as a common drawing object in <GtCanvas> and
   <GtCustomTrack> implementations and supports a variety of drawing operations
   for both text and basic primitive shapes. */
typedef struct GtGraphics GtGraphics;

/* Draws text in black to the right of (<x>,<y>). The coordinate <y> is used as
   a baseline. */
void   gt_graphics_draw_text(GtGraphics*, double x, double y, const char*);
/* Draws text in black to the right of (<x>,<y>). The coordinate <y> is used as
   a baseline. If the text exceeds the margins, it is clipped. */
void   gt_graphics_draw_text_clip(GtGraphics*, double x, double y, const char*);
/* Synonym to <gt_graphics_draw_text()> */
#define gt_graphics_draw_text_left(g,x,y,t) \
        gt_graphics_draw_text(g,x,y,t);
/* Draws text in black centered at (<x>,<y>). The coordinate <y> is used as a
   baseline. */
void   gt_graphics_draw_text_centered(GtGraphics*, double x, double y,
                                      const char*);
/* Draws text in black to the left of (<x>,<y>). The coordinate <y> is used as a
   baseline. */
void   gt_graphics_draw_text_right(GtGraphics*, double x, double y,
                                   const char*);
/* Draws text in a given <GtColor> to the right of (<x>,<y>). The coordinate <y>
   is used as a baseline. */
void   gt_graphics_draw_colored_text(GtGraphics*, double x, double y,
                                     GtColor, const char*);
/* Returns the height of a capital letter in pixels/points. */
double gt_graphics_get_text_height(GtGraphics*);
/* Sets the background color of the <GtGraphics> to a specific color.
   Note that this may only be supported for bitmap output formats. */
int gt_graphics_set_background_color(GtGraphics*, GtColor);
/* Returns the width of the given string in pixels/points. */
double gt_graphics_get_text_width(GtGraphics*, const char *text);
/* Sets basic font family, slant and weight options. Font families are
   implementation-specific, e.g. in Cairo there is no operation to list
   available family names on the system, but the standard CSS2 generic family
   names, ("serif", "sans-serif", "cursive", "fantasy", "monospace"), are
   likely to work as expected.*/
void   gt_graphics_set_font(GtGraphics *g, const char *family,
                            FontSlant slant, FontWeight weight, double size);
/* Returns the width of the image in pixels/points. */
double gt_graphics_get_image_width(GtGraphics*);
/* Returns the height of the image in pixels/points. */
double gt_graphics_get_image_height(GtGraphics*);
/* Set margins (space to the image boundaries that are clear of elements)
   in the graphics.
   <margin_x> denotes the Margin to the left and right, in pixels.
   <margin_y> denotes the Margin to the top and bottom, in pixels. */
void   gt_graphics_set_margins(GtGraphics*, double margin_x,
                                  double margin_y);
/* Returns the horizontal margins in pixels/points. */
double gt_graphics_get_xmargins(GtGraphics*);
/* Returns the vertical margins in pixels/points. */
double gt_graphics_get_ymargins(GtGraphics*);
/* Draws a horizontal line of length <width> beginning at the given coordinates
   to the right in the color <color> with stroke width <stroke_width>. */
void   gt_graphics_draw_horizontal_line(GtGraphics *g, double x, double y,
                                        GtColor color, double width,
                                        double stroke_width);
/* Draws a vertical line of length <length> beginning at the given coordinates
   downwards in the color <color> with stroke width <stroke_width>. */
void   gt_graphics_draw_vertical_line(GtGraphics *g, double x, double y,
                                      GtColor color, double length,
                                      double stroke_width);
/* Draws a line beginning at (<x>,<y>) to (<xto>,<yto>) in the color <color>
   with stroke width <stroke_width>.  */
void   gt_graphics_draw_line(GtGraphics *g, double x, double y,
                             double xto, double yto, GtColor color,
                             double stroke_width);
/* Draws a arrow-like box glyph at (<x>,<y>) where these are the top left
   coordinates. The box extends <width> pixels (incl. arrowhead) into the x
   direction and <height> pixels into the y direction. It will be filled with
   <fill_color> and stroked with width <stroke_width> and color <stroke_color>.
   The width of the arrowhead is given by the <arrow_width> parameter.
   The <arrow_status> parameter determines whether an arrowhead will be drawn
   at the left or right end, both ends, or none.
   If <dashed> is set to true, then the outline will be dashed instead of
   solid.*/
void   gt_graphics_draw_box(GtGraphics*, double x, double y, double width,
                            double height, GtColor fill_color,
                            ArrowStatus arrow_status, double arrow_width,
                            double stroke_width, GtColor stroke_color,
                            bool dashed);
/* Draws a transparent box with a dashed line at the center at (<x>,<y>)
   (where these are the top left coordinates). The box extends <width> pixels
   (incl. arrowhead) into the x direction and <height> pixels into the y
   direction. It will be stroked with width <stroke_width> and color
   <stroke_color>. The width of the arrowhead is given by the <arrow_width>
   parameter. The <arrow_status> parameter determines whether an arrowhead will
   be drawn at the left or right end, both ends, or none. */
void   gt_graphics_draw_dashes(GtGraphics*, double x, double y,
                                  double width, double height,
                                  ArrowStatus arrow_status, double arrow_width,
                                  double stroke_width, GtColor stroke_color);
/* Draws a caret (``hat'') style glyph at (<x>,<y>) (where these are the top
   left coordinates). The box extends <width> pixels (incl. arrowhead) into the
   x direction and <height> pixels into the y direction. It will be stroked
   with width <stroke_width> and color <stroke_color>. The width of the
   arrowhead is given by the <arrow_width> parameter. The <arrow_status>
   parameter determines whether an arrowhead will be drawn at the left or right
   end, both ends, or none. */
void   gt_graphics_draw_caret(GtGraphics*, double x, double y, double width,
                              double height, ArrowStatus arrow_status,
                              double arrow_width,  double stroke_width,
                              GtColor stroke_color);
/* Draws a rectangle at (<x>,<y>) where these are the top left coordinates.
   The rectangle extends <width> pixels (incl. arrowhead) into the x
   direction and <height> pixels into the y direction. It will be filled with
   <fill_color> if <filled> is set to true and stroked with width <stroke_width>
   and color <stroke_color> if <stroked> is set to true. */
void   gt_graphics_draw_rectangle(GtGraphics*, double x, double y,
                                  bool filled, GtColor fill_color,
                                  bool stroked, GtColor stroke_color,
                                  double stroke_width, double width,
                                  double height);
/* Draws an arrowhead at (<x>,<y>) where these are the top left coordinates.
   The direction is determined by the <arrow_status> parameter. */
void   gt_graphics_draw_arrowhead(GtGraphics*, double x, double y, GtColor,
                                  ArrowStatus arrow_status);
/* Draws a curve over the full visible image width (without margins) at
   (<x>,<y>) where these are the top left coordinates. As input, the array of
   double values <data> with <ndata> data points is used. The <valrange> gives
   the minimum and maximum value of the displayed data. If a value outside the
   data range is encountered, the drawing will be stopped at this data point. */
void   gt_graphics_draw_curve_data(GtGraphics *g, double x, double y,
                                   GtColor color,
                                   double data[], GtUword ndata,
                                   GtRange valrange, GtUword height);
/* Write out the <GtGraphics> object to the given file with <filename>. */
int    gt_graphics_save_to_file(const GtGraphics*, const char *filename,
                                GtError*);
/* Write out the <GtGraphics> object to the given <stream>. */
void   gt_graphics_save_to_stream(const GtGraphics*, GtStr *stream);
/* Deletes the the <GtGraphics> object. */
void   gt_graphics_delete(GtGraphics*);

#endif
