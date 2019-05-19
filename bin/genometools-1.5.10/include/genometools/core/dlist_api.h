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

#ifndef DLIST_API_H
#define DLIST_API_H

#include "core/error_api.h"
#include "core/fptr_api.h"
#include "core/types_api.h"

/* A double-linked list which is sorted according to a <GtCompare> compare
   function (<qsort(3)>-like, only if one was supplied to the constructor). */
typedef struct GtDlist GtDlist;
typedef struct GtDlistelem GtDlistelem;

/* Return a new <GtDlist> object sorted according to <compar> function. If
   <compar> equals <NULL>, no sorting is enforced. */
GtDlist*     gt_dlist_new(GtCompare compar);

/* Return a new <GtDlist> object sorted according to <compar> function. If
   <compar> equals <NULL>, no sorting is enforced. Use <data> to supply
   additional data to the comparator function. */
GtDlist*     gt_dlist_new_with_data(GtCompareWithData compar, void *data);

/* Return the first <GtDlistelem> object in <dlist>. */
GtDlistelem* gt_dlist_first(const GtDlist *dlist);

/* Return the last <GtDlistelem> object in <dlist>. */
GtDlistelem* gt_dlist_last(const GtDlist *dlist);

/* Return the first <GtDlistelem> object in <dlist> which contains data
   identical to <data>. Takes O(n) time. */
GtDlistelem* gt_dlist_find(const GtDlist *dlist, void *data);

/* Return the number of <GtDlistelem> objects in <dlist>. */
GtUword      gt_dlist_size(const GtDlist *dlist);

/* Add a new <GtDlistelem> object containing <data> to <dlist>.
   Usually O(n), but O(1) if data is added in sorted order. */
void         gt_dlist_add(GtDlist *dlist, void *data);

/* Remove <dlistelem> from <dlist> and free it. */
void         gt_dlist_remove(GtDlist *dlist, GtDlistelem *dlistelem);

/* Example for usage of the <GtDlist> class. */
int          gt_dlist_example(GtError *err);

/* Delete <dlist>. */
void         gt_dlist_delete(GtDlist *dlist);

/* Return the successor of <dlistelem>, or <NULL> if the element is the last
   one in the <GtDlist>. */
GtDlistelem* gt_dlistelem_next(const GtDlistelem *dlistelem);

/* Return the predecessor of <dlistelem>, or <NULL> if the element is the
   first one in the <GtDlist>. */
GtDlistelem* gt_dlistelem_previous(const GtDlistelem *dlistelem);

/* Return the data pointer attached to <dlistelem>. */
void*        gt_dlistelem_get_data(const GtDlistelem *dlistelem);

#endif
