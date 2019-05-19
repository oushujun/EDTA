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

#ifndef ERROR_API_H
#define ERROR_API_H

#include <stdarg.h>
#include <stdbool.h>
#include "core/assert_api.h"

/*
   This class is used for the handling of ___user errors___ in __GenomeTools__.
   Thereby, the actual <GtError> object is used to store the __error message__
   while it is signaled by the return value of the called function, if an error
   occurred.

   By convention in __GenomeTools__, the <GtError> object is always passed into
   a function as the last parameter and -1 (or <NULL> for constructors) is used
   as return value to indicate that an error occurred.
   Success is usually indicated by 0 as return value or via a non-<NULL> object
   pointer for constructors.

   It is possible to use <NULL> as an <GtError> object, if one is not interested
   in the actual error message.

   Functions which do not get an <GtError> object cannot fail due to a user
   error and it is not necessary to check their return code for an error
   condition.
*/
typedef struct GtError GtError;

/* Return a new <GtError> object */
GtError*   gt_error_new(void);
/* Insert an assertion to check that the error <err> is not set or is <NULL>.
   This macro should be used at the beginning of every routine which has an
   <GtError*> argument to make sure the error propagation has been coded
   correctly. */
#define gt_error_check(err)\
        gt_assert(!err || !gt_error_is_set(err))
/* Set the error message stored in <err> according to <format> (as in
   <printf(3)>).  */
void        gt_error_set(GtError *err, const char *format, ...)
              __attribute__ ((format (printf, 2, 3)));
/* Set the error message stored in <err> according to <format> (as in
   <vprintf(3)>). */
void        gt_error_vset(GtError *err, const char *format, va_list ap);
/* Set the error message stored in <err> to <msg>. */
void        gt_error_set_nonvariadic(GtError *err, const char *msg);
/* Return <true> if the error <err> is set, <false> otherwise. */
bool        gt_error_is_set(const GtError *err);
/* Unset the error <err>. */
void        gt_error_unset(GtError *err);
/* Return the error string stored in <err> (the error must be set). */
const char* gt_error_get(const GtError *err);
/* Delete the error object <err>. */
void        gt_error_delete(GtError *err);

#endif
