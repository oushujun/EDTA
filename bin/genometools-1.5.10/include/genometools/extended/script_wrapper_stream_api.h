/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SCRIPT_WRAPPER_STREAM_API_H
#define SCRIPT_WRAPPER_STREAM_API_H

#include "extended/node_stream_api.h"

/* Implements the <GtScriptWrapperStream> interface. This stream is
   only used to store pointers to external callbacks, e.g. written in a
   scripting language. This class does not store any state or logic, relying
   on the developer of the external custom stream class to do so.  */
typedef struct GtScriptWrapperStream GtScriptWrapperStream;

typedef int (*GtScriptWrapperStreamNextFunc)(GtGenomeNode **gn,
                                             GtError *err);
typedef int (*GtScriptWrapperStreamFreeFunc)(void*);

/* Creates a new <GtScriptWrapperStream> given a next and a free function. */
GtNodeStream* gt_script_wrapper_stream_new(GtScriptWrapperStreamNextFunc,
                                           GtScriptWrapperStreamFreeFunc);

#endif
