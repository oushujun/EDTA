/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef RDB_API_H
#define RDB_API_H

#include "core/error_api.h"
#include "core/phase_api.h"
#include "core/str_api.h"
#include "core/cstr_table_api.h"
#include "core/strand_api.h"

/* the ``GtRDB'' interface for relational database abstraction */
typedef struct GtRDB GtRDB;
typedef struct GtRDBMembers GtRDBMembers;
typedef struct GtRDBStmt GtRDBStmt;

typedef struct GtRDBClass GtRDBClass;
typedef struct GtRDBStmtClass GtRDBStmtClass;

#include "extended/rdb_visitor_api.h"

GtRDB*        gt_rdb_ref(GtRDB*);
GtRDBStmt*    gt_rdb_prepare(GtRDB *db, const char *query,
                             GtUword num_params, GtError *err);
int           gt_rdb_accept(GtRDB *db, GtRDBVisitor *v, GtError *err);
GtCstrTable*  gt_rdb_get_indexes(GtRDB *db, GtError *err);
GtCstrTable*  gt_rdb_get_tables(GtRDB *db, GtError *err);
GtUword gt_rdb_last_inserted_id(GtRDB *db, const char *table, GtError*);
void          gt_rdb_delete(GtRDB*);

int           gt_rdb_stmt_reset(GtRDBStmt *stmt, GtError *err);
int           gt_rdb_stmt_bind_int(GtRDBStmt *stmt, GtUword param_no,
                                   int val, GtError *err);
int           gt_rdb_stmt_bind_ulong(GtRDBStmt *stmt, GtUword param_no,
                                     GtUword val, GtError *err);
int           gt_rdb_stmt_bind_string(GtRDBStmt *stmt, GtUword param_no,
                                      const char *val, GtError *err);
int           gt_rdb_stmt_bind_double(GtRDBStmt *stmt, GtUword param_no,
                                      double val, GtError *err);
int           gt_rdb_stmt_exec(GtRDBStmt *stmt, GtError *err);
int           gt_rdb_stmt_get_ulong(GtRDBStmt *stmt, GtUword field_no,
                                    GtUword *result, GtError *err);
int           gt_rdb_stmt_get_int(GtRDBStmt *stmt, GtUword field_no,
                                  int *result, GtError *err);
int           gt_rdb_stmt_get_string(GtRDBStmt *stmt, GtUword field_no,
                                     GtStr *result, GtError *err);
int           gt_rdb_stmt_get_double(GtRDBStmt *stmt, GtUword field_no,
                                     double *result, GtError *err);
void          gt_rdb_stmt_delete(GtRDBStmt *stmt);
#endif
