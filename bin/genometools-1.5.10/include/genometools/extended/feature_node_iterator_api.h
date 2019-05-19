/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef FEATURE_NODE_ITERATOR_API_H
#define FEATURE_NODE_ITERATOR_API_H

#include "extended/feature_node_api.h"

typedef struct GtFeatureNodeIterator GtFeatureNodeIterator;

/* Return a new <GtFeatureNodeIterator*> which performs a depth-first
   traversal of <feature_node> (including <feature_node> itself).
   It ignores pseudo-features. */
GtFeatureNodeIterator* gt_feature_node_iterator_new(const GtFeatureNode
                                                    *feature_node);
/* Return a new <GtFeatureNodeIterator*> which iterates over all direct
   children of <feature_node> (without <feature_node> itself). */
GtFeatureNodeIterator* gt_feature_node_iterator_new_direct(const GtFeatureNode
                                                           *feature_node);
/* Return the next <GtFeatureNode*> in <feature_node_iterator> or <NULL> if none
   exists. */
GtFeatureNode*         gt_feature_node_iterator_next(GtFeatureNodeIterator
                                                     *feature_node_iterator);
int                    gt_feature_node_iterator_example(GtError*);
/* Delete <feature_node_iterator>. */
void                   gt_feature_node_iterator_delete(GtFeatureNodeIterator
                                                       *feature_node_iterator);

#endif
