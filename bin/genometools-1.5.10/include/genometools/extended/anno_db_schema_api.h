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

#ifndef ANNO_DB_SCHEMA_API_H
#define ANNO_DB_SCHEMA_API_H

#include "extended/feature_index_api.h"
#include "core/error_api.h"
#include "extended/rdb_api.h"

/* The <GtAnnoDBSchema> interface for a database-backed abstract
   <GtFeatureIndex> factory. */
typedef struct GtAnnoDBSchema GtAnnoDBSchema;

/* Returns a <GtFeatureIndex> object representing <GtRDB> object <db>
   interpreted as having schema <schema>. Returns NULL if an error occurred,
   <err> is set accordingly. */
GtFeatureIndex* gt_anno_db_schema_get_feature_index(GtAnnoDBSchema *schema,
                                                    GtRDB *db, GtError *err);

/* Deletes <schema> and frees all associated memory. */
void            gt_anno_db_schema_delete(GtAnnoDBSchema *schema);

#endif
