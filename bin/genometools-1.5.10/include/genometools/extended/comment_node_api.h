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

#ifndef COMMENT_NODE_API_H
#define COMMENT_NODE_API_H

#include "core/error_api.h"

/* Implements the <GtGenomeNode> interface. Comment nodes correspond to comment
   lines in GFF3 files (i.e., lines which start with a single ``<#>''). */
typedef struct GtCommentNode GtCommentNode;

#include "extended/genome_node_api.h"

const GtGenomeNodeClass* gt_comment_node_class(void);

/* Return a new <GtCommentNode> object representing a <comment>. Please note
   that the single leading ``<#>'' which denotes comment lines in GFF3 files
   should not be part of <comment>. */
GtGenomeNode*            gt_comment_node_new(const char *comment);

/* Return the comment stored in <comment_node>. */
const char*              gt_comment_node_get_comment(const GtCommentNode
                                                     *comment_node);

/* Test whether the given genome node is a comment node. If so, a pointer to the
   meta node is returned. If not, NULL is returned. Note that in most cases,
   one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtCommentNode*           gt_comment_node_try_cast(GtGenomeNode *gn);

/* Test whether the given genome node is a comment node. If so, a pointer to the
   meta node is returned. If not, an assertion fails. */
GtCommentNode*           gt_comment_node_cast(GtGenomeNode *gn);

#endif
