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

#ifndef LAYOUT_API_H
#define LAYOUT_API_H

#include "annotationsketch/canvas_api.h"
#include "annotationsketch/diagram_api.h"
#include "annotationsketch/style_api.h"
#include "annotationsketch/text_width_calculator_api.h"
#include "core/range_api.h"

/* The <GtLayout> class represents contents (tracks) of a <GtDiagram> broken up
   into lines such that a given horizontal space allotment given in pixels
   or points is used up most efficiently. This is done using the <GtLineBreaker>
   and <GtTextWidthCalculator> classes. As defaults, Cairo-based instances of
   these classes are used but can be specified separately.

   A <GtLayout> can be queried for the height of the laid out representation and
   finally be rendered to a <GtCanvas>. */
typedef struct GtLayout GtLayout;

/* A function describing the order of tracks based on their track identifier
   strings <s1> and <s2>. Must return a negative value if the track with ID <s1>
   should appear before the track with ID <s2> and a positive value if <s1>
   should appear after <s2>. Returning a value of 0 will result in an undefined
   ordering of <s1> and <s2>. */
typedef int (*GtTrackOrderingFunc)(const char *s1, const char *s2, void *data);

/* Creates a new <GtLayout> object for the contents of <diagram>.
   The layout is done for a target image width of <width> and using the rules in
   <GtStyle> object <style>. */
GtLayout*     gt_layout_new(GtDiagram *diagram, unsigned int width, GtStyle*,
                            GtError*);
/* Like <gt_layout_new()>, but allows use of a different <GtTextWidthCalculator>
   implementation. */
GtLayout*     gt_layout_new_with_twc(GtDiagram*,
                                     unsigned int width,
                                     GtStyle*,
                                     GtTextWidthCalculator*,
                                     GtError*);
/* Sets the <GtTrackOrderingFunc> comparator function <func> which defines an
   order on the tracks contained in <layout>. This determines the order in
   which the tracks are drawn vertically.
   Additional data necessary in the comparator function can be given in <data>,
   the caller is responsible to free it. */
void          gt_layout_set_track_ordering_func(GtLayout *layout,
                                                GtTrackOrderingFunc func,
                                                void *data);
void          gt_layout_unset_track_ordering_func(GtLayout *layout);
/* Calculates the height of <layout> in pixels. The height value is written to
   the location pointed to by <result>. If an error occurs during the
   calculation, this function returns -1 and <err> is set accordingly.
   Returns 0 on success. */
int           gt_layout_get_height(GtLayout *layout,
                                   GtUword *result,
                                   GtError *err);
/* Renders <layout> on the <target_canvas>. */
int           gt_layout_sketch(GtLayout *layout, GtCanvas *target_canvas,
                               GtError*);
/* Destroys a layout. */
void          gt_layout_delete(GtLayout*);

#endif
