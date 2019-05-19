/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef CUSTOM_TRACK_SCRIPT_WRAPPER_API_H
#define CUSTOM_TRACK_SCRIPT_WRAPPER_API_H

#include "annotationsketch/custom_track_api.h"
#include "annotationsketch/graphics_api.h"
#include "annotationsketch/style_api.h"
#include "core/error_api.h"
#include "core/range_api.h"
#include "core/str_api.h"

/* Implements the <GtCustomTrack> interface. This custom track is
   only used to store pointers to external callbacks, e.g. written in a
   scripting language. This class does not store any state, relying on the
   developer of the external custom track class to do so.  */
typedef struct GtCustomTrackScriptWrapper GtCustomTrackScriptWrapper;

typedef int           (*GtCtScriptRenderFunc)(GtGraphics*,
                                              unsigned int,
                                              GtRange*,
                                              GtStyle*,
                                              GtError*);
typedef GtUword (*GtCtScriptGetHeightFunc)(void*);
typedef void          (*GtCtScriptGetTitleFunc)(void*, GtStr*);
typedef void          (*GtCtScriptFreeFunc)(void*);

/* Creates a new <GtCustomTrackScriptWrapper> object. */
GtCustomTrack* gt_custom_track_script_wrapper_new(GtCtScriptRenderFunc
                                                             render_func,
                                                  GtCtScriptGetHeightFunc
                                                             get_height_func,
                                                  GtCtScriptGetTitleFunc
                                                             get_title_func,
                                                  GtCtScriptFreeFunc
                                                             free_func);
#endif
