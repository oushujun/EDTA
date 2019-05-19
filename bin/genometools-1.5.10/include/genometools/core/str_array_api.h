/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STR_ARRAY_API_H
#define STR_ARRAY_API_H

#include "core/str_api.h"

/* <GtStrArray*> objects are arrays of string which grow on demand. */
typedef struct GtStrArray GtStrArray;

/* Return a new <GtStrArray> object. */
GtStrArray*   gt_str_array_new(void);
/* Increases the reference to a GtStrArray. */
GtStrArray*   gt_str_array_ref(GtStrArray*);
/* Add <cstr> to <str_array>. Thereby, an internal copy of <cstr> is created. */
void          gt_str_array_add_cstr(GtStrArray *str_array, const char *cstr);
/* Add the non <\0>-terminated <cstr> with given <length> to <str_array>.
   Thereby, an internal copy of <cstr> is created. */
void          gt_str_array_add_cstr_nt(GtStrArray *str_array, const char *cstr,
                                       GtUword length);
/* Add <str> to <str_array>. Thereby, an internal copy of <str> is created. */
void          gt_str_array_add(GtStrArray *str_array, const GtStr *str);
/* Return pointer to internal string with number <strnum> of <str_array>.
   <strnum> must be smaller than <gt_str_array_size(str_array)>. */
const char*   gt_str_array_get(const GtStrArray *str_array,
                               GtUword strnum);
/* Set the string with number <strnum> in <str_array> to <cstr>. */
void          gt_str_array_set_cstr(GtStrArray *str_array, GtUword strnum,
                                    const char *cstr);
/* Set the string with number <strnum> in <str_array> to <str>. */
void          gt_str_array_set(GtStrArray *str_array, GtUword strnum,
                               const GtStr *str);
/* Set the size of <str_array> to <size>. <size> must be smaller or equal than
   <gt_str_array_size(str_array)>. */
void          gt_str_array_set_size(GtStrArray *str_array, GtUword size);
/* Set the size of <str_array> to 0. */
void          gt_str_array_reset(GtStrArray *str_array);
/* Return the number of strings stored in <str_array>. */
GtUword       gt_str_array_size(const GtStrArray *str_array);
/* Delete <str_array>. */
void          gt_str_array_delete(GtStrArray *str_array);

#endif
