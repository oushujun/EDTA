/*
  Copyright (c) 2010      Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2011-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_ITERATOR_API_H
#define MATCH_ITERATOR_API_H

#include "core/error_api.h"
#include "extended/match_api.h"

typedef struct GtMatchIterator GtMatchIterator;
typedef struct GtMatchIteratorClass GtMatchIteratorClass;

typedef enum {
  GT_MATCHER_STATUS_OK,
  GT_MATCHER_STATUS_END,
  GT_MATCHER_STATUS_ERROR
} GtMatchIteratorStatus;

/* Advances <mp> by one, returning the next match.
   Writes a pointer to the next match to the position pointed to by <match>.
   Returns GT_MATCHER_STATUS_OK when the match could be delivered and there are
   more matches to come, GT_MATCHER_STATUS_END when no more matches are
   available, and GT_MATCHER_STATUS_ERROR if an error occurred. <err> is set
   accordingly. */
GtMatchIteratorStatus gt_match_iterator_next(GtMatchIterator *mp,
                                             GtMatch **match, GtError *err);

/* Deletes <mp>, freeing all associated space. */
void                  gt_match_iterator_delete(GtMatchIterator *mp);

#endif
