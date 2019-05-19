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

#ifndef RDB_SQLITE_API_H
#define RDB_SQLITE_API_H

/* The <GtRDBSqlite> class implements the <GtRDB> interface using the
   SQLite embedded database backend. This implementation is only available if
   compiled with the option ``with-sqlite=yes''. */
typedef struct GtRDBSqlite GtRDBSqlite;

/* The <GtRDBStmtSqlite> class implements the <GtRDBStmt> interface using the
   SQLite embedded database backend. This implementation is only available if
   compiled with the option ``with-sqlite=yes''. */
typedef struct GtRDBStmtSqlite GtRDBStmtSqlite;

#include "extended/rdb_api.h"

/* Creates a new <GtRDBSqlite> object from the database file located at
   <dbpath>. Returns NULL on error, <err> is set accordingly. */
GtRDB* gt_rdb_sqlite_new(const char *dbpath, GtError *err);

#endif
