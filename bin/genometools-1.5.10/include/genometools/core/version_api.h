/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef VERSION_API_H
#define VERSION_API_H

#include <stdarg.h>

/* Version module */

/* Check that the __GenomeTools__ library in use is compatible with the given
   version. Generally you would pass in the constants <GT_MAJOR_VERSION>,
   <GT_MINOR_VERSION>, and <GT_MICRO_VERSION> as the three arguments to this
   function.

   Returns <NULL> if the __GenomeTools__ library is compatible with the given
   version, or a string describing the version mismatch, if the library is not
   compatible.
*/
const char* gt_version_check(unsigned int required_major,
                             unsigned int required_minor,
                             unsigned int required_micro);

/* Return the version of the __GenomeTools__ library in use as a string. */
const char* gt_version(void);

#endif
