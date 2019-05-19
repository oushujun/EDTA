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

#ifndef WARNING_API_H
#define WARNING_API_H

#include <stdarg.h>

/* Warning module */

/* Handler type used to process warnings. */
typedef void (*GtWarningHandler)(void *data, const char *format, va_list ap);

/* Print a warning according to <format> and <...>, if a handler is set. */
void gt_warning(const char *format, ...)
  __attribute__ ((format (printf, 1, 2)));

/* Disable that warnings are shown. That is, subsequent <gt_warning()> calls
   have no effect. */
void gt_warning_disable(void);

/* Set <warn_handler> to handle all warnings issued with <gt_warning()>.
   The <data> is passed to <warning_handler> on each invocation. */
void gt_warning_set_handler(GtWarningHandler warn_handler, void *data);

/* The default warning handler which prints on <stderr>.
   "warning: " is prepended and a newline is appended to the message defined by
   <format> and <ap>. Does not use <data>. */
void gt_warning_default_handler(void *data, const char *format, va_list ap);

/* Return currently used <GtWarningHandler>. */
GtWarningHandler gt_warning_get_handler(void);

/* Return currently used <data> which is passed to the currently used
   <GtWarningHandler>. */
void* gt_warning_get_data(void);

#endif
