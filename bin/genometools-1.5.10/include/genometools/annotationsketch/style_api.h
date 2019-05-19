/*
  Copyright (c) 2007-2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef STYLE_API_H
#define STYLE_API_H

#include <stdbool.h>
#include "annotationsketch/color_api.h"
#include "core/error_api.h"
#include "core/str_api.h"
#include "extended/genome_node_api.h"

typedef enum {
  GT_STYLE_QUERY_OK,
  GT_STYLE_QUERY_NOT_SET,
  GT_STYLE_QUERY_ERROR
} GtStyleQueryStatus;

/* Objects of the <GtStyle> class hold __AnnotationSketch__ style information
   like colors, margins, collapsing options, and others. The class provides
   methods to set values of various types. Each value is organized into
   a __section__ and is identified by a __key__. That is, a __section__, __key__
   pair must uniquely identify a value. */
typedef struct GtStyle GtStyle;

/* Creates a new <GtStyle> object. */
GtStyle*           gt_style_new(GtError*);
/* Increments the reference count of the given <GtStyle>. */
GtStyle*           gt_style_ref(GtStyle*);
/* Enables unsafe mode (``io'' and ``os'' libraries loaded). */
void               gt_style_unsafe_mode(GtStyle*);
/* Enables safe mode (``io'' and ``os'' libraries not accessible). */
void               gt_style_safe_mode(GtStyle*);
/* Returns true if <sty> is in unsafe mode. */
bool               gt_style_is_unsafe(GtStyle *sty);
/* Creates a independent (``deep'') copy of the given <GtStyle> object. */
GtStyle*           gt_style_clone(const GtStyle*, GtError*);
/* Loads and executes Lua style file with given <filename>.
   This file must define a global table called __style__. */
int                gt_style_load_file(GtStyle*, const char *filename, GtError*);
/* Loads and executes Lua style code from the given <GtStr> <instr>.
   This code must define a global table called __style__. */
int                gt_style_load_str(GtStyle*, GtStr *instr, GtError*);
/* Generates Lua code which represents the given <GtStyle> object and
   writes it into the <GtStr> object <outstr>.*/
int                gt_style_to_str(const GtStyle*, GtStr *outstr, GtError*);
/* Reloads the Lua style file. */
void               gt_style_reload(GtStyle*);
/* Sets a color value in the <GtStyle> for section <section> and <key> to a
   certain <color>. */
void               gt_style_set_color(GtStyle*, const char *section,
                                      const char *key, const GtColor *color);
/* Retrieves a color value from <style> for key <key> in section <section>.
   The color is written to the location pointed to by <result>. Optionally, a
   feature node pointer <fn> can be specified for handling in node-specific
   callbacks.
   Because color definitions can be functions, <gt_style_get_color()> can fail
   at runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and
   <err> is set accordingly.
   If the color was not specified in <style>, a grey default color
   is written to <result> and GT_STYLE_QUERY_NOT_SET is returned so the caller
   can provide a custom default.
   In case of successful retrieval of an existing color, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_color(const GtStyle *style, const char *section,
                                      const char *key, GtColor *result,
                                      GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_color(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_color_with_track(const GtStyle *style,
                                                 const char *section,
                                                 const char *key,
                                                 GtColor *result,
                                                 GtFeatureNode *fn,
                                                 const GtStr *track_id,
                                                 GtError *err);
/* Set string with key <key> in <section> to <value>. */
void               gt_style_set_str(GtStyle*, const char *section,
                                    const char *key, GtStr *value);
/* Retrieves a string value from <style> for key <key> in section <section>.
   The string is written to the <GtStr> object <result>, overwriting its prior
   contents. Optionally, a feature node pointer <fn> can be specified for
   handling in node-specific callbacks.
   Because color definitions can be functions, <gt_style_get_str()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the string was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing string, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_str(const GtStyle *style, const char *section,
                                    const char *key, GtStr *result,
                                    GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_str(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_str_with_track(const GtStyle *style,
                                               const char *section,
                                               const char *key,
                                               GtStr *result,
                                               GtFeatureNode *fn,
                                               const GtStr *track_id,
                                               GtError *err);
/* Set numeric value of key <key> in <section> to <number>. */
void               gt_style_set_num(GtStyle*, const char *section,
                                    const char *key, double number);
/* Retrieves a numeric value from <style> for key <key> in section <section>.
   The value is written to the location pointed to by <result>. Optionally, a
   feature node pointer <fn> can be specified for handling in node-specific
   callbacks.
   Because the definitions can be functions, <gt_style_get_num()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the number was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing number, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_num(const GtStyle *style, const char *section,
                                    const char *key, double *result,
                                    GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_num(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_num_with_track(const GtStyle *style,
                                               const char *section,
                                               const char *key,
                                               double *result,
                                               GtFeatureNode *fn,
                                               const GtStr *track_id,
                                               GtError *err);
/* Set boolean value of key <key> in <section> to <val>. */
void               gt_style_set_bool(GtStyle*, const char *section,
                                     const char *key, bool val);
/* Retrieves a boolean value from <style> for key <key> in section <section>.
   The value is written to the location pointed to by <result>. Optionally, a
   feature node pointer <fn> can be specified for handling in node-specific
   callbacks.
   Because the definitions can be functions, <gt_style_get_bool()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the value was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing boolean, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_bool(const GtStyle *style, const char *section,
                                     const char *key, bool *result,
                                     GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_bool(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_bool_with_track(const GtStyle *style,
                                                const char *section,
                                                const char *key,
                                                bool *result,
                                                GtFeatureNode *fn,
                                                const GtStr *track_id,
                                                GtError *err);
/* Unset value of key <key> in <section>. */
void               gt_style_unset(GtStyle*, const char *section,
                                  const char *key);
/* Deletes this <style>. */
void               gt_style_delete(GtStyle *style);

#endif
