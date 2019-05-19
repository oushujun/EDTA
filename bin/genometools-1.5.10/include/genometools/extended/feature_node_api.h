/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_NODE_API_H
#define FEATURE_NODE_API_H

#include "core/phase_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/str_array_api.h"

/* Implements the <GtGenomeNode> interface. A single feature node corresponds
   to a GFF3 feature line (i.e., a line which does not start with <#>).
   Part-of relationships (which are realized in GFF3 with the <Parent> and <ID>
   attributes) are realized in the C API with the
   <gt_feature_node_add_child()> method.

   Besides the ``mere'' feature nodes two ``special'' feature nodes exist:
   multi-features and pseudo-features.

   Multi-features represent features which span multiple lines (it is indicated
   in GFF3 files by the fact, that each line has the same ID attribute).

   To check if a feature is a multi-feature use the method
   <gt_feature_node_is_multi()>.
   Multi-features are connected via a ``representative''. That is, two features
   are part of the same multi-feature if they have the same representative.
   The feature node representative can be be retrieved via the
   <gt_feature_node_get_multi_representative()> method.

   Pseudo-features became a technical necessity to be able to pass related
   top-level features as a single entity through the streaming machinery.
   There are two cases in which a pseudo-feature has to be introduced.

   First, if a multi-feature has no parent. In this case all features which
   comprise the multi-feature become the children of a pseudo-feature.

   Second, if two or more top-level features have the same children (and are
   thereby connected). In this case all these top-level features become the
   children of a pseudo-feature.

   It should be clear from the explanation above that pseudo-features make only
   sense as top-level features (a fact which is enforced in the code).

   Pseudo-features are typically ignored during a traversal to give the illusion
   that they do not exist. */
typedef struct GtFeatureNode GtFeatureNode;

#include "extended/genome_node_api.h"

/* Return an new <GtFeatureNode> object on sequence with ID <seqid> and type
   <type> which lies from <start> to <end> on strand <strand>.
   The <GtFeatureNode*> stores a new reference to <seqid>, so make sure you do
   not modify the original <seqid> afterwards!
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode*  gt_feature_node_new(GtStr *seqid, const char *type,
                                   GtUword start, GtUword end,
                                   GtStrand strand);

/* Return a new pseudo-<GtFeatureNode> object on sequence with ID <seqid> which
   lies from <start> to <end> on strand <strand>. Pseudo-features do not have a
   type.  The <GtFeatureNode > stores a new reference to <seqid>, so make sure
   you do not modify the original <seqid> afterwards.  <start> and <end> always
   refer to the forward strand, therefore <start> has to be smaller or equal
   than <end>. */
GtGenomeNode*  gt_feature_node_new_pseudo(GtStr *seqid, GtUword start,
                                          GtUword end, GtStrand strand);

/* Return a new pseudo-<GtFeatureNode> object which uses <feature_node> as
   template.  That is, the sequence ID, range, strand, and source are taken from
   <feature_node>. */
GtGenomeNode*  gt_feature_node_new_pseudo_template(GtFeatureNode *feature_node);

/* Return the ``standard gene'' (mainly for testing purposes). */
GtGenomeNode*  gt_feature_node_new_standard_gene(void);

/* Add <child> feature node to <parent> feature node.
   <parent> takes ownership of <child>.*/
void           gt_feature_node_add_child(GtFeatureNode *parent,
                                         GtFeatureNode *child);

/* Return the source of <feature_node>. If no source has been set, "." is
   returned. Corresponds to column 2 of GFF3 feature lines. */
const char*    gt_feature_node_get_source(const GtFeatureNode *feature_node);

/* Set the <source> of <feature_node>. Stores a new reference to <source>.
   Corresponds to column 2 of GFF3 feature lines. */
void           gt_feature_node_set_source(GtFeatureNode *feature_node,
                                          GtStr *source);

/* Return <true> if <feature_node> has a defined source (i.e., on different
   from "."). <false> otherwise. */
bool           gt_feature_node_has_source(const GtFeatureNode *feature_node);

/* Return the type of <feature_node>. Corresponds to column 3 of GFF3 feature
   lines. */
const char*    gt_feature_node_get_type(const GtFeatureNode *feature_node);

/* Set the type of <feature_node> to <type>. */
void           gt_feature_node_set_type(GtFeatureNode *feature_node,
                                        const char *type);

/* Return <true> if <feature_node> has given <type>, <false> otherwise. */
bool           gt_feature_node_has_type(GtFeatureNode *feature_node,
                                        const char *type);

/* Return the number of children for given <feature_node>. */
GtUword        gt_feature_node_number_of_children(const GtFeatureNode
                                                  *feature_node);

/* Return the number of children of type <node>
   for given GtFeatureNode <parent>. */
GtUword        gt_feature_node_number_of_children_of_type(const GtFeatureNode
                                                          *parent,
                                                          const GtFeatureNode
                                                          *node);

/* Return <true> if the score of <feature_node> is defined, <false>
   otherwise. */
bool           gt_feature_node_score_is_defined(const GtFeatureNode
                                                *feature_node);
/* Return the score of <feature_node>. The score has to be defined.
   Corresponds to column 6 of GFF3 feature lines. */
float          gt_feature_node_get_score(const GtFeatureNode *feature_node);

/* Set the score of <feature_node> to <score>. */
void           gt_feature_node_set_score(GtFeatureNode *feature_node,
                                         float score);

/* Unset the score of <feature_node>. */
void           gt_feature_node_unset_score(GtFeatureNode *feature_node);

/* Return the strand of <feature_node>. Corresponds to column 7 of GFF3
   feature lines. */
GtStrand       gt_feature_node_get_strand(const GtFeatureNode *feature_node);

/* Set the strand of <feature_node> to <strand>. */
void           gt_feature_node_set_strand(GtFeatureNode *feature_node,
                                          GtStrand strand);

/* Return the phase of <feature_node>. Corresponds to column 8 of GFF3 feature
   lines. */
GtPhase        gt_feature_node_get_phase(const GtFeatureNode *feature_node);

/* Set the phase of <feature_node> to <phase>. */
void           gt_feature_node_set_phase(GtFeatureNode *feature_node,
                                         GtPhase phase);

/* Return the attribute of <feature_node> with the given <name>.
   If no such attribute has been added, <NULL> is returned.
   The attributes are stored in column 9 of GFF3 feature lines. */
const char*    gt_feature_node_get_attribute(const GtFeatureNode *feature_node,
                                             const char *name);

/* Return a string array containing the used attribute names of <feature_node>.
   The caller is responsible to free the returned <GtStrArray*>. */
GtStrArray*    gt_feature_node_get_attribute_list(const GtFeatureNode
                                                  *feature_node);

/* Add attribute <tag>=<value> to <feature_node>. <tag> and <value> must at
   least have length 1. <feature_node> must not contain an attribute with the
   given <tag> already. You should not add Parent and ID attributes, use
   <gt_feature_node_add_child()> to denote part-of relationships. */
void           gt_feature_node_add_attribute(GtFeatureNode *feature_node,
                                             const char *tag,
                                             const char *value);

/* Set attribute <tag> to new <value> in <feature_node>, if it exists already.
   Otherwise the attribute <tag>=<value> is added to <feature_node>.
   <tag> and <value> must at least have length 1.
   You should not set Parent and ID attributes, use
   <gt_feature_node_add_child()> to denote part-of relationships. */
void           gt_feature_node_set_attribute(GtFeatureNode* feature_node,
                                             const char *tag,
                                             const char *value);

/* Remove attribute <tag> from <feature_node>. <feature_node> must contain an
   attribute with the given <tag> already! You should not remove Parent and ID
   attributes. */
void           gt_feature_node_remove_attribute(GtFeatureNode* feature_node,
                                                const char *tag);
/* Delivers the key (in <attr_name>)and value (in <attr_value>) of an
   attribute. The <data> parameter carries over arbitrary user data from the
   <gt_feature_node_foreach_attribute> call. */

typedef void (*GtFeatureNodeAttributeIterFunc)(const char *attr_name,
                                               const char *attr_value,
                                               void *data);

/* Calls <func> for each attribute in <feature_node>. Use <data> to forward
   arbitrary data during traversal. */
void           gt_feature_node_foreach_attribute(GtFeatureNode *feature_node,
                                            GtFeatureNodeAttributeIterFunc func,
                                            void *data);

/* Return <true> if <feature_node> is a multi-feature, <false> otherwise. */
bool           gt_feature_node_is_multi(const GtFeatureNode *feature_node);

/* Return <true> if <feature_node> is a pseudo-feature, <false> otherwise. */
bool           gt_feature_node_is_pseudo(const GtFeatureNode *feature_node);

/* Make <feature_node> the representative of a multi-feature.
   Thereby <feature_node> becomes a multi-feature. */
void           gt_feature_node_make_multi_representative(GtFeatureNode
                                                         *feature_node);

/* Set the multi-feature representative of <feature_node> to <representative>.
   Thereby <feature_node> becomes a multi-feature. */
void           gt_feature_node_set_multi_representative(GtFeatureNode
                                                        *feature_node,
                                                        GtFeatureNode
                                                        *representative);

/* Unset the multi-feature status of <feature_node> and remove its multi-feature
   representative. */
void           gt_feature_node_unset_multi(GtFeatureNode *feature_node);

/* Return the representative of the multi-feature <feature_node>. */
GtFeatureNode* gt_feature_node_get_multi_representative(GtFeatureNode
                                                        *feature_node);

/* Returns <true>, if the given <feature_node_a> has the same seqid, feature
   type, range, strand, and phase as <feature_node_b>.
   Returns <false> otherwise. */
bool           gt_feature_node_is_similar(const GtFeatureNode *feature_node_a,
                                          const GtFeatureNode *feature_node_b);

/* Marks the given <feature_node>. */
void           gt_feature_node_mark(GtFeatureNode*);

/* If the given <feature_node> is marked it will be unmarked. */
void           gt_feature_node_unmark(GtFeatureNode*);

/* Returns <true> if the given <feature_node> graph contains a marked node. */
bool           gt_feature_node_contains_marked(GtFeatureNode *feature_node);

/* Returns <true> if the (top-level) <feature_node> is marked. */
bool           gt_feature_node_is_marked(const GtFeatureNode *feature_node);

/* Removes the parent-child relationship between the leaf node <leafn> and a
   parent node <tree>. <tree> needs not be the direct parent of the <leafn>.
   Note that <leafn> is freed, use <gt_genome_node_ref()> to increase the
   reference count if deletion is not wanted. As a side effect, the tree
   property of <tree> is set back to an undefined state. */
void           gt_feature_node_remove_leaf(GtFeatureNode *tree,
                                           GtFeatureNode *leafn);

/* Test whether the given genome node is a feature node. If so, a pointer to the
   feature node is returned. If not, NULL is returned. Note that in most cases,
   one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtFeatureNode* gt_feature_node_try_cast(GtGenomeNode *gn);

/* Test whether the given genome node is a feature node. If so, a pointer to the
   feature node is returned. If not, an assertion fails. */
GtFeatureNode* gt_feature_node_cast(GtGenomeNode *gn);

#endif
