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

#ifndef SELECT_STREAM_API_H
#define SELECT_STREAM_API_H

#include "core/strand_api.h"
#include "extended/node_stream_api.h"

/* Implements the <GtNodeStream> interface. A <GtSelectStream> selects certain
   nodes it retrieves from its node source and passes them along. */
typedef struct GtSelectStream GtSelectStream;

typedef int (*GtSelectNodeFunc)(GtGenomeNode *gn, void *data, GtError *err);

/* Create a <GtSelectStream> object which selects genome nodes it retrieves from
   its <in_stream> and passes them along if they meet the criteria defined by
   the other arguments. All comment nodes are selected.
   If <seqid> is defined, a genome node must have it to be selected.
   If <source> is defined, a genome node must have it to be selected.
   If <contain_range> is defined, a genome node must be contained in it to be
   selected.
   If <overlap_range> is defined, a genome node must overlap it to be selected.
   If <strand> is defined, a (top-level) genome node must have it to be
   selected.
   If <targetstrand> is defined, a feature with a target attribute must have
   exactly one of it and its strand must equal <targetstrand>.
   If <had_cds> is <true>, all top-level features are selected which have a
   child with type __CDS__.
   If <max_gene_length> is defined, only genes up to the this length are
   selected.
   If <max_gene_num> is defined, only so many genes are selected.
   If <min_gene_score> is defined, only genes with at least this score are
   selected.
   If <max_gene_score> is defined, only genes with at most this score are
   selected.
   If <min_average_splice_site_prob> is defined, feature nodes which have
   splice sites must have at least this average splice site score to be
   selected.
   If <feature_num> is defined, just the <feature_num>th feature node occurring
   in the <in_stream> is selected.
   If <select_files> is defined and has at least one entry, the entries are
   evaluated as Lua scripts containing functions taking <GtGenomeNodes> that
   are evaluated to boolean values to determine selection. <select_logic>
   can be "OR" or "AND", defining how the results from the select scripts are
   combined.
   Returns a pointer to a new <GtSelectStream> or NULL on error (<err> is set
   accordingly).
*/
GtNodeStream* gt_select_stream_new(GtNodeStream *in_stream,
                                   GtStr *seqid,
                                   GtStr *source,
                                   const GtRange *contain_range,
                                   const GtRange *overlap_range,
                                   GtStrand strand,
                                   GtStrand targetstrand,
                                   bool has_CDS,
                                   GtUword max_gene_length,
                                   GtUword max_gene_num,
                                   double min_gene_score,
                                   double max_gene_score,
                                   double min_average_splice_site_prob,
                                   GtUword feature_num,
                                   GtStrArray *select_files,
                                   GtStr *select_logic,
                                   GtError *err);

/* Sets <fp> as a handler function to be called for every <GtGenomeNode> not
   selected by <sstr>. The void pointer <data> can be used for arbitrary user
   data.
*/
void          gt_select_stream_set_drophandler(GtSelectStream *sstr,
                                               GtSelectNodeFunc fp,
                                               void *data);

#endif
