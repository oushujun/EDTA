/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef BLOCK_API_H
#define BLOCK_API_H

#include "core/range_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "extended/feature_node_api.h"

/* The <GtBlock> class represents a portion of screen space which relates to
   a specific ``top-level'' feature (and maybe its collapsed child features).
   It is the smallest layoutable unit in AnnotationSketch and has a caption
   (which may be displayed above the block rendering). */
typedef struct GtBlock GtBlock;

/* Creates a new <GtBlock> object. */
GtBlock*              gt_block_new(void);
/* Increases the reference count. */
GtBlock*              gt_block_ref(GtBlock*);
/* Create a new GtBlock object, setting block parameters (such as strand,
   range) from a given <node> template. */
GtBlock*              gt_block_new_from_node(GtFeatureNode *node);
/* Returns the base range of the <GtBlock>'s top level element. */
GtRange               gt_block_get_range(const GtBlock*);
/* Returns a pointer to the base range of the <GtBlock>'s top level element. */
GtRange*              gt_block_get_range_ptr(const GtBlock *block);
/* Checks whether a <GtBlock> is occupied completely by a single element. */
bool                  gt_block_has_only_one_fullsize_element(const GtBlock*);
/* Merges the contents of two <GtBlock>s into the first one. */
void                  gt_block_merge(GtBlock*, GtBlock*);
/* Returns an independent copy of a <GtBlock>. */
GtBlock*              gt_block_clone(GtBlock*);
/* Set whether a block caption should be displayed or not. */
void                  gt_block_set_caption_visibility(GtBlock*, bool);
/* Returns whether a block caption should be displayed or not. */
bool                  gt_block_caption_is_visible(const GtBlock*);
/* Sets the <GtBlock>'s caption to <caption>. */
void                  gt_block_set_caption(GtBlock*, GtStr *caption);
/* Returns the <GtBlock>'s caption. */
GtStr*                gt_block_get_caption(const GtBlock*);
/* Sets the <GtBlock>'s strand to <strand>. */
void                  gt_block_set_strand(GtBlock*, GtStrand strand);
/* Returns the <GtBlock>'s strand. */
GtStrand              gt_block_get_strand(const GtBlock*);
/* Returns the <GtBlock>'s top level feature as a <GtFeatureNode> object. */
GtFeatureNode*        gt_block_get_top_level_feature(const GtBlock*);
/* Returns the number of elements in the <GtBlock>. */
GtUword         gt_block_get_size(const GtBlock*);
/* Returns the feature type of the <GtBlock>. */
const char*           gt_block_get_type(const GtBlock*);
/* Deletes a <GtBlock>. */
void                  gt_block_delete(GtBlock*);

#endif
