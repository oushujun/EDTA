/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef SCRIPT_FILTER_API_H
#define SCRIPT_FILTER_API_H

#include "core/error_api.h"
#include "extended/feature_node_api.h"

typedef struct GtScriptFilter GtScriptFilter;

GtScriptFilter* gt_script_filter_new(const char *file, GtError *err);
GtScriptFilter* gt_script_filter_new_unsafe(const char *file, GtError *err);
GtScriptFilter* gt_script_filter_new_from_string(const char *script_string,
                                                 GtError *err);
GtScriptFilter* gt_script_filter_ref(GtScriptFilter *script_filter);
const char*     gt_script_filter_get_name(GtScriptFilter *script_filter,
                                          GtError *err);
const char*     gt_script_filter_get_description(GtScriptFilter *script_filter,
                                                 GtError *err);
const char*     gt_script_filter_get_short_description(
                                                  GtScriptFilter *script_filter,
                                                  GtError *err);
const char*     gt_script_filter_get_author(GtScriptFilter *script_filter,
                                            GtError *err);
const char*     gt_script_filter_get_version(GtScriptFilter *script_filter,
                                             GtError *err);
const char*     gt_script_filter_get_email(GtScriptFilter *script_filter,
                                           GtError *err);
bool            gt_script_filter_validate(GtScriptFilter *script_filer,
                                          GtError *err);
int             gt_script_filter_run(GtScriptFilter *script_filter,
                                     GtFeatureNode *gf,
                                     bool *select_node,
                                     GtError *err);
void            gt_script_filter_delete(GtScriptFilter *script_filter);

#endif
