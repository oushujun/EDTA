/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef REGION_NODE_API_H
#define REGION_NODE_API_H

/* Implements the <GtGenomeNode> interface. Region nodes correspond to the
   <##sequence-region> lines in GFF3 files.*/
typedef struct GtRegionNode GtRegionNode;

#include "extended/genome_node_api.h"
#include "core/str_api.h"

const GtGenomeNodeClass* gt_region_node_class(void);

/* Create a new <GtRegionNode*> representing sequence with ID <seqid> from
   base position <start> to base position <end> (1-based).
   <start> has to be smaller or equal than <end>.
   The <GtRegionNode*> stores a new reference to <seqid>, so make sure you do
   not modify the original <seqid> afterwards! */
GtGenomeNode* gt_region_node_new(GtStr *seqid, GtUword start,
                                               GtUword end);

/* Test whether the given genome node is a region node. If so, a pointer to the
   region node is returned. If not, NULL is returned. Note that in most cases,
   one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtRegionNode* gt_region_node_try_cast(GtGenomeNode *gn);

/* Test whether the given genome node is a region node. If so, a pointer to the
   region node is returned. If not, an assertion fails. */
GtRegionNode* gt_region_node_cast(GtGenomeNode *gn);

#endif
