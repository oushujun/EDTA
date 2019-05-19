/*
  Copyright (c) 2012, 2015 Gordon Gremme <gordon@gremme.org>

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

#ifndef META_NODE_API_H
#define META_NODE_API_H

#include "core/error_api.h"

/* Implements the <GtGenomeNode> interface. Meta nodes correspond to meta
   lines in GFF3 files (i.e., lines which start with  ``<##>'') which are not
   sequence-region lines. */
typedef struct GtMetaNode GtMetaNode;

#include "extended/genome_node_api.h"

const GtGenomeNodeClass* gt_meta_node_class(void);

/* Return a new <GtMetaNode> object representing a <meta_directive> with the
   corresponding <meta_data>. Please note that the leading ``<##>'' which
   denotes meta lines in GFF3 files should not be part of the
   <meta_directive>.
   The <meta_directive> must not be <NULL>, the <meta_data> can be <NULL>. */
GtGenomeNode*            gt_meta_node_new(const char *meta_directive,
                                          const char *meta_data);

/* Return the meta directive stored in <meta_node>. */
const char*              gt_meta_node_get_directive(const GtMetaNode
                                                    *meta_node);

/* Return the meta data stored in <meta_node>. Can return NULL! */
const char*              gt_meta_node_get_data(const GtMetaNode *meta_node);

/* Test whether the given genome node is a meta node. If so, a pointer to the
   meta node is returned. If not, NULL is returned. Note that in most cases,
   one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtMetaNode*              gt_meta_node_try_cast(GtGenomeNode *gn);

/* Test whether the given genome node is a meta node. If so, a pointer to the
   meta node is returned. If not, an assertion fails. */
GtMetaNode*              gt_meta_node_cast(GtGenomeNode *gn);

#endif
