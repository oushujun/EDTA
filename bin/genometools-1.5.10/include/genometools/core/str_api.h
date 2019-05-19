/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef STR_API_H
#define STR_API_H

#include <stdio.h>
#include "core/error_api.h"
#include "core/types_api.h"

/* Objects of the <GtStr> class are strings which grow on demand. */
typedef struct GtStr GtStr;

/* Return an empty <GtStr> object. */
GtStr*  gt_str_new(void);
/* Return a new <GtStr> object whose content is set to <cstr>. */
GtStr*  gt_str_new_cstr(const char *cstr);
/* Return a clone of <str>. */
GtStr*  gt_str_clone(const GtStr *str);
/* Increase the reference count for <str> and return it.
   If <str> is <NULL>, <NULL> is returned without any side effects. */
GtStr*  gt_str_ref(GtStr *str);
/* Return the content of <str>.  Never returns NULL, and the content is always
   <\0>-terminated */
char*   gt_str_get(const GtStr *str);
/* Set the content of <str> to <cstr>. */
void    gt_str_set(GtStr *str, const char *cstr);
/* Append the string <src> to <dest>. */
void    gt_str_append_str(GtStr *dest, const GtStr *src);
/* Append the <\0>-terminated <cstr> to <str>. */
void    gt_str_append_cstr(GtStr *str, const char *cstr);
/* Append the (not necessarily <\0>-terminated) <cstr> with given <length> to
   <str>. */
void    gt_str_append_cstr_nt(GtStr *str, const char *cstr, GtUword length);
/* Append character <c> to <str>. */
void    gt_str_append_char(GtStr *str, char c);
/* Append double <d> to <str> with given <precision>. */
void    gt_str_append_double(GtStr *str, double d, int precision);
/* Append double <d> to <str> in scientific notation e.g. 0.52e10,
   with given <precision>. */
void gt_str_append_sci_double(GtStr *dest, double d, int precision);
/* Append <ulong> to <str>. */
GT_DEPRECATED("use gt_str_append_uword() instead")
void    gt_str_append_ulong(GtStr *str, GtUword ulong);
/* Append <ulong> to <str>. */
void    gt_str_append_uword(GtStr *str, GtUword uword);
/* Append <intval> to <str>. */
void    gt_str_append_int(GtStr *str, int intval);
/* Append <uint> to <str>. */
void    gt_str_append_uint(GtStr *str, unsigned int uint);
/* Set length of <str> to <length>. <length> must be smaller or equal than
   <gt_str_length(str)>. */
void    gt_str_set_length(GtStr *str, GtUword length);
/* Reset <str> to length 0. */
void    gt_str_reset(GtStr *str);
/* Compare <str1> and <str2> and return the result (similar to <strcmp(3)>). */
int     gt_str_cmp(const GtStr *str1, const GtStr *str2);
/* Return the length of <str>. If <str> is <NULL>, 0 is returned. */
GtUword gt_str_length(const GtStr *str);
/* Decrease the reference count for <str> or delete it, if this was the last
   reference. */
void    gt_str_delete(GtStr *str);

#endif
