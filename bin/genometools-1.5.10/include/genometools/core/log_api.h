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

#ifndef LOG_API_H
#define LOG_API_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

/* Log module */

/* Enable logging. */
void gt_log_enable(void);

/* Returns true if logging is enabled, false otherwise */
bool gt_log_enabled(void);

/* Prints the log message obtained from format and following parameters
   according if logging is enabled. The logging output is prefixed with the
   string "debug: " and finished by a newline.  */
void  gt_log_log(const char *format, ...)
  __attribute__ ((format (printf, 1, 2)));

/* Prints the log message obtained from format and following parameter according
   to if logging is enabled analog to gt_log_log(). But in contrast to
   gt_log_log() gt_log_vlog() does not accept individual arguments but a single
   va_list argument instead. */
void  gt_log_vlog(const char *format, va_list);

/* Return logging file pointer. */
FILE* gt_log_fp(void);

/* Set logging file pointer to <fp>. */
void  gt_log_set_fp(FILE *fp);

#endif
