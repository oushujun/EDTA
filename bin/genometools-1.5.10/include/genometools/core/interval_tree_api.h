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

#ifndef INTERVAL_TREE_API_H
#define INTERVAL_TREE_API_H

#include "core/array_api.h"
#include "core/fptr_api.h"

/* This is an interval tree data structure, implemented according to
   Cormen et al., Introduction to Algorithms, 2nd edition, MIT Press,
   Cambridge, MA, USA, 2001 */
typedef struct GtIntervalTree GtIntervalTree;
typedef struct GtIntervalTreeNode GtIntervalTreeNode;

typedef int (*GtIntervalTreeIteratorFunc)(GtIntervalTreeNode*, void*);

/* Creates a new <GtIntervalTreeNode>. Transfers ownership of <data> to interval
   tree if inserted into a <GtIntervalTree> in which a
   <GtIntervalTreeDataFreeFunc> is set. */
GtIntervalTreeNode* gt_interval_tree_node_new(void *data,
                                              GtUword low,
                                              GtUword high);

/* Returns a pointer to the data associated with node <node>. */
void*               gt_interval_tree_node_get_data(GtIntervalTreeNode* node);

/* Creates a new <GtIntervalTree>. If a <GtFree> function is given as an
   argument, it is applied on the data pointers in all inserted nodes when the
   <GtIntervalTree> is deleted. */
GtIntervalTree*     gt_interval_tree_new(GtFree);

/* Returns the number of elements in the <GtIntervalTree>. */
GtUword       gt_interval_tree_size(GtIntervalTree*);

/* Returns the first node in the <GtIntervalTree> which overlaps the given
   range (from <start> to <end>). */
GtIntervalTreeNode* gt_interval_tree_find_first_overlapping(GtIntervalTree*,
                                                            GtUword start,
                                                            GtUword end);

/* Inserts node <node> into <tree>. */
void                gt_interval_tree_insert(GtIntervalTree *tree,
                                            GtIntervalTreeNode *node);

/* Collects data pointers of all <GtIntervalTreeNode>s in the tree which
   overlap with the query range (from <start> to <end>) in a <GtArray>. */
void                gt_interval_tree_find_all_overlapping(GtIntervalTree*,
                                                          GtUword start,
                                                          GtUword end,
                                                          GtArray*);

/* Call <func> for all <GtIntervalTreeNode>s in the tree which overlap with
   the query range (from <start> to <end>). Use <data> to pass in arbitrary
   user data. */
void                gt_interval_tree_iterate_overlapping(GtIntervalTree *it,
                                                GtIntervalTreeIteratorFunc func,
                                                GtUword start,
                                                GtUword end,
                                                void *data);

/* Traverses the <GtIntervalTree> in a depth-first fashion, applying <func> to
   each node encountered. The <data> pointer can be used to reference arbitrary
   data needed in the <GtIntervalTreeIteratorFunc>. */
int                 gt_interval_tree_traverse(GtIntervalTree*,
                                              GtIntervalTreeIteratorFunc func,
                                              void *data);

/* Removes the entry referenced by <node> from the <GtIntervalTree>.
   The data attached to <node> is freed according to the free function defined
   in the tree.
   Note that the memory pointed to by <node> can be re-used internally,
   referencing other data in the tree. Make sure to handle this pointer as
   expired after calling <gt_interval_tree_remove()>! */
void                gt_interval_tree_remove(GtIntervalTree*,
                                            GtIntervalTreeNode *node);

/* Deletes a <GtIntervalTree>. If a <GtFree> function was set in the tree
   constructor, data pointers specified in the nodes are freed using the given
   <GtFree> function. */
void                gt_interval_tree_delete(GtIntervalTree*);

#endif
