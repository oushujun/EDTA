/*
  Copyright (c) 2014 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef DEPRECATED_API_H
#define DEPRECATED_API_H

/* Deprecated module */

/* Deprecated functions, typedefs and structs in API headers should be annotated
   with this macro to get warnings when the gcc or clang compiler is used.
   <msg> should inform about an alternative API. */
#if defined(__clang__)
#if __has_extension(attribute_deprecated_with_message)
#define GT_DEPRECATED(msg) \
         __attribute__((deprecated (msg)))
#else
#define GT_DEPRECATED(msg) \
         __attribute__((deprecated))
#endif
/* not clang with messages for deprecated */
#elif defined(__GNUC__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 405)
#define GT_DEPRECATED(msg) \
         __attribute__((deprecated (msg)))
/* not gcc >= 4.5 */
#elif defined(__GNUC__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 301)
#define GT_DEPRECATED(msg) \
         __attribute__((deprecated))
#else /* not gcc >= 3.1 */
#define GT_DEPRECATED(msg)
#endif

#endif
