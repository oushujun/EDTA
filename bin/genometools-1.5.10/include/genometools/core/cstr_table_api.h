/*
  Copyright (c) 2008-2010, 2012 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008            Center for Bioinformatics, University of Hamburg

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

#ifndef CSTR_TABLE_API_H
#define CSTR_TABLE_API_H

#include "core/error_api.h"
#include "core/str_array_api.h"

/* Implements a table of C strings. */
typedef struct GtCstrTable GtCstrTable;

/* Return a new <GtCstrTable> object. */
GtCstrTable*  gt_cstr_table_new(void);
/* Add <cstr> to <table>. <table> must not already contain <cstr>! */
void          gt_cstr_table_add(GtCstrTable *table, const char *cstr);
/* If a C string equal to <cstr> is contained in <table>, it is returned.
   Otherwise <NULL> is returned. */
const char*   gt_cstr_table_get(const GtCstrTable *table, const char *cstr);
/* Return a <GtStrArray*> which contains all <cstr>s added to <table> in
   alphabetical order. The caller is responsible to free it! */
GtStrArray*   gt_cstr_table_get_all(const GtCstrTable *table);
/* Remove <cstr> from <table>. */
void          gt_cstr_table_remove(GtCstrTable *table, const char *cstr);
/* Reset <table> (that is, remove all contained C strings). */
void          gt_cstr_table_reset(GtCstrTable *table);
/* Delete C string <table>. */
void          gt_cstr_table_delete(GtCstrTable *table);

#endif
