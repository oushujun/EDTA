/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef ANNO_DB_GFFLIKE_API_H
#define ANNO_DB_GFFLIKE_API_H

/* The <GtAnnoDBGFFlike> class implements the <GtAnnoDBSchema> interface,
   using a database schema specifically tailored to store __GenomeTools__
   annotations. */
typedef struct GtAnnoDBGFFlike GtAnnoDBGFFlike;

#include "core/error_api.h"
#include "extended/anno_db_schema_api.h"

/* Creates a new <GtAnnoDBGFFlike> schema object. */
GtAnnoDBSchema* gt_anno_db_gfflike_new(void);

/* Retrieves all features contained in <gfi> into <results>. Returns 0 on
   success, a negative value otherwise. The message in <err> is set
   accordingly. */
int             gt_feature_index_gfflike_get_all_features(GtFeatureIndex *gfi,
                                                          GtArray *results,
                                                          GtError *err);

int             gt_anno_db_gfflike_unit_test(GtError *err);

#endif
