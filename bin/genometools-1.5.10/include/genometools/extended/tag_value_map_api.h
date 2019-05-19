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

#ifndef TAG_VALUE_MAP_API_H
#define TAG_VALUE_MAP_API_H

#include "core/error_api.h"

/* A very simple tag/value map absolutely optimized for space (i.e., memory
   consumption) on the cost of time. Basically, each read/write access costs
   O(n) time, whereas n denotes the accumulated length of all tags and values
   contained in the map. Tags and values cannot have length 0.

   The implementation as a char* shines through (also to save one additional
   memory allocation), therefore the usage is a little bit different compared
   to other __GenomeTools__ classes.
   See the implementation of <gt_tag_value_map_example()> for an ussage
   example. */
typedef char* GtTagValueMap;

/* Iterator function used to iterate over tag/value maps. A <tag>/<value> pair
   and user <data> are given as arguments. */
typedef void (*GtTagValueMapIteratorFunc)(const char *tag, const char *value,
                                          void *data);

/* Return a new <GtTagValueMap> object which stores the given <tag>/<value>
   pair. */
GtTagValueMap gt_tag_value_map_new(const char *tag, const char *value);
/* Add <tag>/<value> pair to <tag_value_map>. <tag_value_map> must not contain
   the given <tag> already! */
void          gt_tag_value_map_add(GtTagValueMap *tag_value_map,
                                   const char *tag,
                                   const char *value);
/* Set the given <tag> in <tag_value_map> to <value>. */
void          gt_tag_value_map_set(GtTagValueMap *tag_value_map,
                                   const char *tag, const char *value);
/* Return value corresponding to <tag> from <tag_value_map>. If <tag_value_map>
   does not contain such a value, <NULL> is returned. */
const char*   gt_tag_value_map_get(const GtTagValueMap tag_value_map,
                                   const char *tag);

/* Return the number of tag-value pairs in <tag_value_map>. */
GtUword       gt_tag_value_map_size(const GtTagValueMap tag_value_map);

/* Removes the given <tag> from <tag_value_map>. <tag_value_map> must contain
   the given <tag> already! Also, at least one tag-value pair must remain in
   the map. */
void          gt_tag_value_map_remove(GtTagValueMap *tag_value_map,
                                      const char *tag);

/* Apply <iterator_func> to each tag/value pair contained in <tag_value_map> and
   pass <data> along. */
void          gt_tag_value_map_foreach(const GtTagValueMap tag_value_map,
                                       GtTagValueMapIteratorFunc iterator_func,
                                       void *data);
/* Implements an example useage of a tag/value map. */
int           gt_tag_value_map_example(GtError *err);
/* Delete <tag_value_map>. */
void          gt_tag_value_map_delete(GtTagValueMap tag_value_map);

#endif
