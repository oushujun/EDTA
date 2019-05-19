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

#ifndef REC_MAP_API_H
#define REC_MAP_API_H

#include "extended/feature_node_api.h"

/* A <GtRecMap> object contains a mapping from a 2D coordinate pair
   which identifies a rectangle in a rendered image to the <GtFeatureNode> it
  represents. The rectangle is defined by the coordinates of its upper left
  (``northwest'') and lower right (``southeast'') points.

  <GtRecMap> objects are created by an <GtImageInfo> object which is filled
  during the generation of an image by __AnnotationSketch__. */
typedef struct GtRecMap GtRecMap;

/* Creates a new <GtRecMap> for feature <f> with the given coordinates. */
GtRecMap*            gt_rec_map_new(double nw_x, double nw_y, double se_x,
                                    double se_y, GtFeatureNode *f);
/* Increases the reference count of <rm>. */
GtRecMap*            gt_rec_map_ref(GtRecMap *rm);
/* Retrieve __x__ value of the the upper left point of the rectangle. */
double               gt_rec_map_get_northwest_x(const GtRecMap*);
/* Retrieve __y__ value of the the upper left point of the rectangle. */
double               gt_rec_map_get_northwest_y(const GtRecMap*);
/* Retrieve __x__ value of the the lower right point of the rectangle. */
double               gt_rec_map_get_southeast_x(const GtRecMap*);
/* Retrieve __y__ value of the the lower right point of the rectangle. */
double               gt_rec_map_get_southeast_y(const GtRecMap*);
/* Retrieve <GtFeatureNode> associated with this rectangle. */
const GtFeatureNode* gt_rec_map_get_genome_feature(const GtRecMap*);
/* Returns <true> if the rectangle represents a block root whose elements
   have not been drawn due to size restrictions. */
bool                 gt_rec_map_has_omitted_children(const GtRecMap*);
/* Deletes a <GtRecMap> and frees all associated memory. */
void                 gt_rec_map_delete(GtRecMap*);

#endif
