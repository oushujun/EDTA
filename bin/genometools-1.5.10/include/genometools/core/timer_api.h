/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef TIMER_API_H
#define TIMER_API_H

#include <stdarg.h>
#include <stdio.h>
#include <stdarg.h>
#include "core/str_api.h"

/* The <GtTimer> class encapsulates a timer which can be used for run-time
   measurements. */
typedef struct GtTimer GtTimer;

/* Return a new <GtTimer> object. */
GtTimer* gt_timer_new(void);
/* Return a new <GtTimer> object with the first <description>. */
GtTimer* gt_timer_new_with_progress_description(const char* description);
/* Start the time measurement on <timer>. */
void     gt_timer_start(GtTimer *timer);
/* Stop the time measurement on <timer>. */
void     gt_timer_stop(GtTimer *timer);
/* Output the current state of <timer> in the format
   ""GT_WD".%06lds real "GT_WD"s user "GT_WD"s system" to file
   pointer <fp> (see <gt_timer_show_formatted>).
   The timer is then stopped. */
void     gt_timer_show(GtTimer *timer, FILE *fp);
/* Output the current state of <timer> in a user-defined format given by <fmt>.
   <fmt> must be a format string for four "GT_WD" numbers, which are filled
   with: elapsed seconds, elapsed microseconds, used usertime in seconds, system
   time in seconds. The output is written to <fp>. The timer is then stopped. */
void     gt_timer_show_formatted(GtTimer *timer, const char *fmt, FILE *fp);
/* return usec of time from start to stop of giben timer. The timer is then
   stopped. */
GtWord gt_timer_elapsed_usec(GtTimer *t);
/* Like <gt_timer_show_formatted()>, but appends the output to <str>. */
void     gt_timer_get_formatted(GtTimer *t, const char *fmt, GtStr *str);
/* Output the current state of <timer> on <fp> since the last call of
   <gt_timer_show_progress()> or the last start of <timer>, along with the
   current description. The timer is not stopped, but updated with <desc> to be
   the next description. */
void     gt_timer_show_progress(GtTimer *timer, const char *desc, FILE *fp);
/* Like <gt_timer_show_progress()>, but allows one to format the description in
   a <printf()>-like fashion. */
void     gt_timer_show_progress_formatted(GtTimer *timer, FILE *fp,
                                          const char *desc, ...);
/* Like <gt_timer_show_progress()>, but allows one to format the description in
   a <vprintf()>-like fashion using a va_list argument <ap>. */
void     gt_timer_show_progress_va(GtTimer *timer, FILE *fp, const char *desc,
                                   va_list ap);
/* Output the overall time measured with <timer> from start to now on <fp>. */
void     gt_timer_show_progress_final(GtTimer *timer, FILE *fp);
/* Show also user and sys time in output of
   <gt_timer_show_progress[_final]()>. */
void     gt_timer_show_cpu_time_by_progress(GtTimer *timer);
/* Hide output of last stage time in <gt_timer_show_progress_final()>. */
void     gt_timer_omit_last_stage(GtTimer *timer);
/* Delete <timer>. */
void     gt_timer_delete(GtTimer *timer);

#endif
