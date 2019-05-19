/*
  Copyright (c) 2007-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007      Malte Mader <mader@zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_INDEX_API_H
#define FEATURE_INDEX_API_H

#include "core/range_api.h"
#include "core/str_array_api.h"
#include "extended/feature_node_api.h"
#include "extended/region_node_api.h"

/* This interface represents a searchable container for <GtFeatureNode> objects,
   typically root nodes of larger structures. How storage and searching takes
   place is left to the discretion of the implementing class.

   Output from a <gt_feature_index_get_features_*()> method should always
   be sorted by feature start position. */
typedef struct GtFeatureIndex GtFeatureIndex;

/* Add <region_node> to <feature_index>. */
int         gt_feature_index_add_region_node(GtFeatureIndex *feature_index,
                                             GtRegionNode *region_node,
                                             GtError *err);
/* Add <feature_node> to <feature_index>, associating it with a sequence
   region denoted by its identifier string. */
int         gt_feature_index_add_feature_node(GtFeatureIndex *feature_index,
                                              GtFeatureNode *feature_node,
                                              GtError *err);
/* Removes node <genome_node> from <feature_index>. */
int         gt_feature_index_remove_node(GtFeatureIndex *feature_index,
                                         GtFeatureNode *node,
                                         GtError *err);
/* Add all features contained in <gff3file> to <feature_index>, if <gff3file> is
   valid. Otherwise, <feature_index> is not changed and <err> is set. */
int         gt_feature_index_add_gff3file(GtFeatureIndex *feature_index,
                                          const char *gff3file, GtError *err);
/* Returns an array of <GtFeatureNodes> associated with a given sequence region
   identifier <seqid>. */
GtArray*    gt_feature_index_get_features_for_seqid(GtFeatureIndex*,
                                                    const char *seqid,
                                                    GtError *err);
/* Look up genome features in <feature_index> for sequence region <seqid> in
   <range> and store them in <results>. */
int         gt_feature_index_get_features_for_range(GtFeatureIndex
                                                    *feature_index,
                                                    GtArray *results,
                                                    const char *seqid,
                                                    const GtRange *range,
                                                    GtError*);
/* Returns the first sequence region identifier added to <feature_index>. */
char*       gt_feature_index_get_first_seqid(const GtFeatureIndex
                                             *feature_index,
                                             GtError *err);
/* Returns a <GtStrArray> of all sequence region identifiers contained in
   <feature_index> (in alphabetical order). */
GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex *feature_index,
                                        GtError *err);
/* Writes the range of all features contained in the <feature_index> for
   region identifier <seqid> to the <GtRange> pointer <range>. */
int         gt_feature_index_get_range_for_seqid(GtFeatureIndex *feature_index,
                                                 GtRange *range,
                                                 const char *seqid,
                                                 GtError *err);
/* Writes the range of the whole sequence region contained in the
   <feature_index> for region identifier <seqid> to the <GtRange>
   pointer <range>. */
int         gt_feature_index_get_orig_range_for_seqid(GtFeatureIndex
                                                                *feature_index,
                                                      GtRange *range,
                                                      const char *seqid,
                                                      GtError *err);
/* Returns <has_seqid>  to true if the sequence region identified by <seqid>
   has been registered in the <feature_index>. */
int         gt_feature_index_has_seqid(const GtFeatureIndex *feature_index,
                                       bool *has_seqid,
                                       const char *seqid,
                                       GtError *err);
/* Calls the save function for the given <feature_index>. The save functions
   must be defined, otherwise the method fails with an assertion (for example,
   the memory based feature index does not have a save function) . */
int         gt_feature_index_save(GtFeatureIndex *feature_index, GtError *err);
/* Deletes the <feature_index> and all its referenced features. */
void        gt_feature_index_delete(GtFeatureIndex*);

#endif
