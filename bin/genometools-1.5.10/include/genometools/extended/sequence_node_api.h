/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef SEQUENCE_NODE_API_H
#define SEQUENCE_NODE_API_H

/* Implements the <GtGenomeNode> interface. Sequence nodes correspond to
   singular embedded FASTA sequences in GFF3 files. */
typedef struct GtSequenceNode GtSequenceNode;

#include "core/str_api.h"
#include "extended/genome_node_api.h"

const GtGenomeNodeClass* gt_sequence_node_class(void);

/* Create a new <GtSequenceNode*> representing a FASTA entry with the given
   <description> and <sequence>. Takes ownership of <sequence>. */
GtGenomeNode*            gt_sequence_node_new(const char *description,
                                              GtStr *sequence);
/* Return the description of <sequence_node>. */
const char*              gt_sequence_node_get_description(const
                                                          GtSequenceNode
                                                          *sequence_node);
/* Return the sequence of <sequence_node>. */
const char*              gt_sequence_node_get_sequence(const GtSequenceNode
                                                       *sequence_node);
/* Return the sequence length of <sequence_node>. */
GtUword                  gt_sequence_node_get_sequence_length(const
                                                              GtSequenceNode
                                                              *sequence_node);

/* Test whether the given genome node is a sequence node. If so, a pointer to
   the sequence node is returned. If not, NULL is returned. Note that in most
   cases, one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtSequenceNode*          gt_sequence_node_try_cast(GtGenomeNode *gn);

/* Test whether the given genome node is a sequence node. If so, a pointer to
   the sequence node is returned. If not, an assertion fails. */
GtSequenceNode*          gt_sequence_node_cast(GtGenomeNode *gn);

#endif
