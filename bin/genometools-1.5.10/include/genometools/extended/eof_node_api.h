/*
  Copyright (c) 2010 Gordon Gremme <gordon@gremme.org>

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

#ifndef EOF_NODE_API_H
#define EOF_NODE_API_H

#include "core/error_api.h"

/* Implements the <GtGenomeNode> interface. EOF nodes mark the barrier between
   separate input files in an GFF3 stream. */
typedef struct GtEOFNode GtEOFNode;

#include "extended/genome_node_api.h"

const GtGenomeNodeClass* gt_eof_node_class(void);

/* Create a new <GtEOFNode*> representing an EOF marker. */
GtGenomeNode*            gt_eof_node_new(void);

/* Test whether the given genome node is an EOF node. If so, a pointer to the
   EOF node is returned. If not, NULL is returned. Note that in most cases,
   one should implement a GtNodeVisitor to handle processing of different
   GtGenomeNode types. */
GtEOFNode*               gt_eof_node_try_cast(GtGenomeNode *gn);

#endif
