/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef FPTR_API_H
#define FPTR_API_H

/* FunctionPointer module */

/* Functions of this type return less than 0 if <a> is __smaller__ than <b>,
   0 if <a> is __equal__ to <b>, and greater 0 if <a> is __larger__ than <b>.
   Thereby, the operators __smaller__, __equal__, and __larger__ are
   implementation dependent.
   Do not count on these functions to return -1, 0, or 1!  */
typedef int  (*GtCompare)(const void *a, const void *b);
/* Similar to <GtCompare>, but with an additional <data> pointer. */
typedef int  (*GtCompareWithData)(const void*, const void*, void *data);
/* The generic free function pointer type. */
typedef void (*GtFree)(void*);

#endif
