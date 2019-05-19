/*
  Copyright (c) 2007-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef MA_API_H
#define MA_API_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* MemoryAllocation module */

/* Allocate ___uninitialized___ space for an object whose size is specified by
   <size> and return it.
   Besides the fact that it never returns <NULL> analog to <malloc(3)>. */
#define       gt_malloc(size)\
              gt_malloc_mem(size, __FILE__, __LINE__)
void*         gt_malloc_mem(size_t size, const char *src_file, int src_line);
/* Allocate contiguous space for an array of <nmemb> objects, each of whose size
   is <size>.  The space is initialized to zero.
   Besides the fact that it never returns <NULL> analog to <calloc(3)>. */
#define       gt_calloc(nmemb, size)\
              gt_calloc_mem(nmemb, size, __FILE__, __LINE__)
void*         gt_calloc_mem(size_t nmemb, size_t size,
                            const char *src_file, int src_line);
/* Change the size of the object pointed to by <ptr> to <size> bytes and return
   a pointer to the (possibly moved) object.
   Besides the fact that it never returns <NULL> analog to <realloc(3)>. */
#define       gt_realloc(ptr, size)\
              gt_realloc_mem(ptr, size, __FILE__, __LINE__)
void*         gt_realloc_mem(void *ptr, size_t size,
                             const char *src_file, int src_line);
/* Free the space pointed to by <ptr>. If <ptr> equals <NULL>, no action occurs.
   Analog to <free(3)>. */
#define       gt_free(ptr)\
              gt_free_mem(ptr, __FILE__, __LINE__)
void          gt_free_mem(void *ptr, const char *src_file, int src_line);
/* Analog to <gt_free()>, but usable as a function pointer. */
void          gt_free_func(void *ptr);

#endif
