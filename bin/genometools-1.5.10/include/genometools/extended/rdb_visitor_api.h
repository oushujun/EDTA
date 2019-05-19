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

#ifndef RDB_VISITOR_API_H
#define RDB_VISITOR_API_H

/* The <GtRDBVisitor> interface, a visitor for <GtRDB> objects. */
typedef struct GtRDBVisitor GtRDBVisitor;

#include "extended/rdb_mysql_api.h"
#include "extended/rdb_sqlite_api.h"

/* Visit a SQLite database <rdbs> with <rdbv>. Returns 0 on success,
   a negative value otherwise, and <err> is set accordingly. */
int   gt_rdb_visitor_visit_sqlite(GtRDBVisitor *rdbv, GtRDBSqlite *rdbs,
                                  GtError *err);
/* Visit a MySQL database <rdbm> with <rdbv>. Returns 0 on success,
   a negative value otherwise, and <err> is set accordingly. */
int   gt_rdb_visitor_visit_mysql(GtRDBVisitor *rdbv, GtRDBMySQL *rdbm,
                                 GtError *err);

/* Delete <rdbv>. */
void  gt_rdb_visitor_delete(GtRDBVisitor *rdbv);

#endif
