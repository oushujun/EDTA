/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef BASENAME_API_H
#define BASENAME_API_H

#include "core/error_api.h"

/* POSIX module */

/*
  This module implements the function gt_basename() according to the
  specifications in
  http://www.unix-systems.org/onlinepubs/7908799/xsh/basename.html
  and
  http://www.opengroup.org/onlinepubs/009695399/

  gt_basename() is equivalent to the function basename(3) which
  is available on most unix systems, but in different libraries and
  with slightly different functionality.

  gt_basename() takes the pathname pointed to by <path> and returns a pointer to
  the final component of the pathname, deleting any trailing '/' characters.

  If <path> consists entirely of the '/' character,  then gt_basename() returns
  a pointer to the string "/".

  On Windows '\' is used instead of '/'.

  If <path> is a null pointer or points to an empty string, gt_basename()
  returns a pointer to the string ".".

  See the implementation of gt_basename_unit_test() for additional examples.

  The caller is responsible for freeing the received pointer!
*/
char* gt_basename(const char *path);
int   gt_basename_unit_test(GtError*);

#endif
