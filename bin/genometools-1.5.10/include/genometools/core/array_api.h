/*
  Copyright (c) 2005-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ARRAY_API_H
#define ARRAY_API_H

#include <stdlib.h>
#include "core/fptr_api.h"
#include "core/types_api.h"

/* <GtArray> objects are generic arrays for elements of a certain size which
   grow on demand. */
typedef struct GtArray GtArray;

/* Return a new <GtArray> object whose elements have the size <size_of_elem>. */
GtArray* gt_array_new(size_t size_of_elem);
/* Increase the reference count for <array> and return it.
   If <array> is <NULL>, <NULL> is returned without any side effects. */
GtArray* gt_array_ref(GtArray *array);
/* Return a clone of <array>. */
GtArray* gt_array_clone(const GtArray *array);
/* Return pointer to element number <index> of <array>. <index> has to be
   smaller than <gt_array_size(array)>. */
void*    gt_array_get(const GtArray *array, GtUword index);
/* Return pointer to first element of <array>. */
void*    gt_array_get_first(const GtArray *array);
/* Return pointer to last element of <array>. */
void*    gt_array_get_last(const GtArray *array);
/* Return pointer to last element of <array> and remove it from <array>. */
void*    gt_array_pop(GtArray *array);
/* Return pointer to the internal space of <array> where the elements are
   stored.  */
void*    gt_array_get_space(const GtArray *array);
/* Add element <elem> to <array>. The size of <elem> must equal the given
   element size when the <array> was created and is determined automatically
   with the <sizeof> operator. */
#define  gt_array_add(array, elem) \
         gt_array_add_elem(array, &(elem), sizeof (elem))
/* Add element <elem> with size <size_of_elem> to <array>. <size_of_elem> must
   equal the given element size when the <array> was created. Usually, this
   method is not used directly and the macro <gt_array_add()> is used
   instead. */
void     gt_array_add_elem(GtArray *array, void *elem, size_t size_of_elem);
/* Add all elements of array <src> to the array <dest>. The element sizes of
   both arrays must be equal. */
void     gt_array_add_array(GtArray *dest, const GtArray *src);
/* Remove element with number <index> from <array> in O(<gt_array_size(array)>)
   time. <index> has to be smaller than <gt_array_size(array)>. */
void     gt_array_rem(GtArray *array, GtUword index);
/* Remove elements starting with number <frompos> up to (and including) <topos>
   from <array> in O(<gt_array_size(array)>) time. <frompos> has to be smaller
   or equal than <topos> and both have to be smaller than
   <gt_array_size(array)>. */
void     gt_array_rem_span(GtArray *array, GtUword frompos,
                                 GtUword topos);
/* Reverse the order of the elements in <array>. */
void     gt_array_reverse(GtArray *array);
/* Set the size of <array> to <size>. <size> must be smaller or equal than
   <gt_array_size(array)>. */
void     gt_array_set_size(GtArray *array, GtUword size);
/* Reset the <array>. That is, afterwards the array has size 0. */
void     gt_array_reset(GtArray *array);
/* Return the size of the elements stored in <array>. */
size_t   gt_array_elem_size(const GtArray *array);
/* Return the number of elements in <array>. If <array> equals <NULL>, 0 is
   returned. */
GtUword  gt_array_size(const GtArray *array);
/* Sort <array> with the given compare function <compar>. */
void     gt_array_sort(GtArray *array, GtCompare compar);
/* Sort <array> in a stable way with the given compare function <compar>. */
void     gt_array_sort_stable(GtArray *array, GtCompare compar);
/* Sort <array> with the given compare function <compar>. Passes a pointer with
   userdata <data> to <compar>. */
void     gt_array_sort_with_data(GtArray *array, GtCompareWithData compar,
                                 void *data);
/* Sort <array> in a stable way with the given compare function <compar>. Passes
   a pointer with userdata <data> to <compar>. */
void     gt_array_sort_stable_with_data(GtArray *array,
                                        GtCompareWithData compar, void *data);
/* Compare the content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same <gt_array_size()> and
   <gt_array_elem_size()>. */
int      gt_array_cmp(const GtArray *array_a, const GtArray *array_b);
/* Decrease the reference count for <array> or delete it, if this was the last
   reference. */
void     gt_array_delete(GtArray *array);

#endif
