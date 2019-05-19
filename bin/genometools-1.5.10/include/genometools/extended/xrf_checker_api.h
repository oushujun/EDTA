/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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

#ifndef XRF_CHECKER_API_H
#define XRF_CHECKER_API_H

#include "core/error_api.h"

/* The <GtXRFChecker> interface, allows one to check the validity of
   Dbxref and Ontology_type attributes. */
typedef struct GtXRFChecker GtXRFChecker;

/* Create a new <GtXRFChecker> from the definitions found in <file_path>.
   Returns NULL on error, and <err> is set accordingly. */
GtXRFChecker* gt_xrf_checker_new(const char *file_path, GtError *err);
/* Increase reference count for <xrf_checker> */
GtXRFChecker* gt_xrf_checker_ref(GtXRFChecker *xrf_checker);
/* Return <true> if <value> is valid for the given <xrf_checker>,
   <false> otherwise. In case of <value> being invalid, <err> is
   set accordingly. */
bool          gt_xrf_checker_is_valid(GtXRFChecker *xrf_checker,
                                      const char *value, GtError *err);
/* Decrease the reference count for <xrf_checker> or delete it. */
void          gt_xrf_checker_delete(GtXRFChecker *xrf_checker);

#endif
