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

#ifndef IMAGE_INFO_API_H
#define IMAGE_INFO_API_H

#include "annotationsketch/rec_map_api.h"
#include "core/error_api.h"

/* The <GtImageInfo> class is a container for 2D coordinate to <GtFeatureNode>
   mappings which could, for example, be used to associate sections of a
   rendered image with GUI widgets or HTML imagemap areas. This information is
   given in the form of <GtRecMap> objects. They are created during the
   image rendering process and stored inside a <GtImageInfo> object for later
   retrieval. Additionally, the rendered width of an image can be obtained via
   a <GtImageInfo> method. */
typedef struct GtImageInfo GtImageInfo;

/* Creates a new <GtImageInfo> object. */
GtImageInfo*     gt_image_info_new(void);
/* Returns the height of the rendered image (in pixels or points). */
unsigned int     gt_image_info_get_height(GtImageInfo *image_info);
/* Returns the total number of mappings in <image_info>. */
GtUword    gt_image_info_num_of_rec_maps(GtImageInfo *image_info);
/* Returns the <i>-th <GtRecMap> mapping in <image_info>. */
const GtRecMap*  gt_image_info_get_rec_map(GtImageInfo *image_info,
                                           GtUword i);
/* Deletes <image_info> and all the <GtRecMap> objects created by it. */
void             gt_image_info_delete(GtImageInfo *image_info);

#endif
