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

#ifndef RDB_MYSQL_API_H
#define RDB_MYSQL_API_H

/* The <GtRDBMySQL> class implements the <GtRDB> interface using the
   MySQL client backend. This implementation is only available if
   compiled with the option ``with-mysql=yes''. */
typedef struct GtRDBMySQL GtRDBMySQL;

/* The <GtRDBStmtMySQL> class implements the <GtRDBStmt> interface using the
   MySQL client database backend. This implementation is only available if
   compiled with the option ``with-mysql=yes''. */
typedef struct GtRDBStmtMySQL GtRDBStmtMySQL;

#include "extended/rdb_api.h"

/* Creates a new <GtRDBSqlite> object from the database accessible on <server>
   port <port>, selecting the database <database> using the credentials given
   by <username> and <password>. Returns NULL on error, <err> is set
   accordingly. */
GtRDB* gt_rdb_mysql_new(const char *server, unsigned int port,
                        const char *database, const char *username,
                        const char *password, GtError *err);

#endif
