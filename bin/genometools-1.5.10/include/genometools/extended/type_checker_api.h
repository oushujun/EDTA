/*
  Copyright (c) 2008, 2011-2012 Gordon Gremme <gordon@gremme.org>
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

#ifndef TYPE_CHECKER_API_H
#define TYPE_CHECKER_API_H

/* The <GtTypeChecker> interface, allows one to check the validity of (genome
   feature) types. */
typedef struct GtTypeChecker GtTypeChecker;

/* Increase the reference count for <type_checker> and return it. */
GtTypeChecker* gt_type_checker_ref(GtTypeChecker *type_checker);
/* Return description of <type_checker>. */
const char*    gt_type_checker_description(GtTypeChecker *type_checker);
/* Return <true> if <type> is a valid type for the given <type_checker>, <false>
   otherwise. */
bool           gt_type_checker_is_valid(GtTypeChecker *type_checker,
                                        const char *type);
/* Return <true> if <child_type> is partof <parent_type>, <false> otherwise. */
bool           gt_type_checker_is_partof(GtTypeChecker *type_checker,
                                         const char *parent_type,
                                         const char *child_type);
/* Return <true> if <child_type> is a <parent_type>, <false> otherwise. */
bool           gt_type_checker_is_a(GtTypeChecker *type_checker,
                                    const char *parent_type,
                                    const char *child_type);
/* Decrease the reference count for <type_checker> or delete it, if this was the
   last reference. */
void           gt_type_checker_delete(GtTypeChecker *type_checker);

#endif
