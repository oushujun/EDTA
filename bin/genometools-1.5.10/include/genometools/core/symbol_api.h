/*
  Copyright (c) 2008-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#ifndef SYMBOL_API_H
#define SYMBOL_API_H

/* Symbol module */

/* Return a symbol (a canonical representation) for <cstr>. An advantage of
   symbols is that they can be compared for equality by a simple pointer
   comparison, rather than using <strcmp()> (as it is done in <gt_strcmp()>).
   Furthermore, a symbol is stored only once in memory for equal <cstr>s, but
   keep in mind that this memory can never be freed safely during the lifetime
   of the calling program. Therefore, it should only be used for a small set of
   <cstr>s. */
const char* gt_symbol(const char *cstr);

#endif
