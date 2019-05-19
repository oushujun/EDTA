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

#ifndef FEATURE_INDEX_MEMORY_API_H
#define FEATURE_INDEX_MEMORY_API_H

#include "extended/feature_index_api.h"

/* The <GtFeatureIndexMemory> class implements a <GtFeatureIndex> in memory.
   Features are organized by region node. Each region node collects its
   feature nodes in an interval tree structure, which allows for efficient
   range queries. */
typedef struct GtFeatureIndexMemory GtFeatureIndexMemory;

/* Creates a new <GtFeatureIndexMemory> object. */
GtFeatureIndex* gt_feature_index_memory_new(void);

/* Returns <ptr> if it is a valid node indexed in <GtFeatureIndexMemory>.
   Otherwise NULL is returned and <err> is set accordingly. */
GtFeatureNode*  gt_feature_index_memory_get_node_by_ptr(GtFeatureIndexMemory*,
                                                        GtFeatureNode *ptr,
                                                        GtError *err);

#endif
