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

#ifndef PHASE_API_H
#define PHASE_API_H

enum GtPhase {
  GT_PHASE_ZERO,     /* '0' */
  GT_PHASE_ONE,      /* '1' */
  GT_PHASE_TWO,      /* '2' */
  GT_PHASE_UNDEFINED /* '.' */
};

/* This enum type defines the possible phases. The following phases are
   defined: <GT_PHASE_ZERO>, <GT_PHASE_ONE>, <GT_PHASE_TWO>, and
   <GT_PHASE_UNDEFINED>. */
typedef enum GtPhase GtPhase;

/* Use this string to map phase enum types to their corresponding character. */
#define GT_PHASE_CHARS \
        "012."

/* Map <phase_char> to the corresponding phase enum type.
   An assertion will fail if <phase_char> is not a valid one. */
GtPhase gt_phase_get(char phase_char);

#endif
