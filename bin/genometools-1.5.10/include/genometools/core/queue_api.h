/*
  Copyright (c) 2008, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2008       Center for Bioinformatics, University of Hamburg

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

#ifndef QUEUE_API_H
#define QUEUE_API_H

#include <stdio.h>
#include "core/error_api.h"
#include "core/types_api.h"

/* <GtQueue> objects are generic queues which can be used to process objects of
   any type in an First-In-First-Out (FIFO) fashion. */
typedef struct GtQueue GtQueue;

/* Return a new <GtQueue> object. */
GtQueue*      gt_queue_new(void);
/* Add <elem> to <queue> (__enqueue__ in computer science terminology). */
void          gt_queue_add(GtQueue *queue, void *elem);
/* Remove the first element from non-empty <queue> and return it (__dequeue__ in
   computer science terminology). */
void*         gt_queue_get(GtQueue *queue);
/* Return the first element in non-empty <queue> without removing it. */
void*         gt_queue_head(GtQueue *queue);
/* Remove <elem> from <queue> (<elem> has to be in <queue>).
   Thereby <queue> is traversed in reverse order, leading to
   O(<gt_queue_size(queue)>) worst-case running time. */
void          gt_queue_remove(GtQueue *queue, void *elem);
/* Return the number of elements in <queue>. */
GtUword gt_queue_size(const GtQueue *queue);
/* Delete <queue>. Elements contained in <queue> are not freed! */
void          gt_queue_delete(GtQueue *queue);

#endif
