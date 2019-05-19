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

#ifndef UNDEF_API_H
#define UNDEF_API_H

#include <float.h>
#include <limits.h>

/* Undef module */

/* The undefined <bool> value. */
#define GT_UNDEF_BOOL \
        (bool) ~0

/* The undefined <char> value. */
#define GT_UNDEF_CHAR \
        CHAR_MAX

/* The undefined <double> value. */
#define GT_UNDEF_DOUBLE \
        DBL_MAX

/* The undefined <float> value. */
#define GT_UNDEF_FLOAT \
        FLT_MAX

/* The undefined <int> value. */
#define GT_UNDEF_INT \
        INT_MIN

/* The undefined <GtWord> value. */
#define GT_UNDEF_WORD \
        LONG_MIN

/* The undefined <long> value. deprecated */
#define GT_UNDEF_LONG \
        GT_UNDEF_WORD

/* The undefined <unsigned char> value. */
#define GT_UNDEF_UCHAR \
        UCHAR_MAX

/* The undefined <unsigned int> value. */
#define GT_UNDEF_UINT \
        ~0U

/* The undefined <GtUword> value. */
#define GT_UNDEF_UWORD \
        ~0UL

/* The undefined <unsigned long> value. deprecated */
#define GT_UNDEF_ULONG \
        GT_UNDEF_UWORD

#endif
