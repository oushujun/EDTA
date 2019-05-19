/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef LOGGER_API_H
#define LOGGER_API_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct GtLogger GtLogger;

/* Creates a new <GtLogger>, with logging <enabled> or not,
   and prefixing all log entries with <prefix> (e.g. "debug").
   The log output is terminated by a newline. All log output will
   be written to <target>. */
GtLogger* gt_logger_new(bool enabled, const char *prefix, FILE *target);
/* Enable logging on <logger>. */
void      gt_logger_enable(GtLogger *logger);
/* Disable logging on <logger>. */
void      gt_logger_disable(GtLogger *logger);
/* Return <true> if logging is enabled on <logger>, false otherwise. */
bool      gt_logger_enabled(GtLogger *logger);
/* Return logging target of <logger>. */
FILE*     gt_logger_target(GtLogger *logger);
/* Set logging target of <logger> to <fp>. */
void      gt_logger_set_target(GtLogger *logger, FILE *fp);
/* Log to target regardless of logging status. */
void      gt_logger_log_force(GtLogger *logger, const char *format, ...)
          __attribute__ ((format (printf, 2, 3)));
/* Log to target depending on logging status. */
void      gt_logger_log(GtLogger *logger, const char *format, ...)
          __attribute__ ((format (printf, 2, 3)));
/* Log to target regardless of logging status, using a va_list argument. */
void      gt_logger_log_va_force(GtLogger *logger, const char *format, va_list);
/* Log to target depending on logging status, using a va_list argument. */
void      gt_logger_log_va(GtLogger *logger, const char *format, va_list);
/* Delete <logger>. */
void      gt_logger_delete(GtLogger *logger);

#endif
