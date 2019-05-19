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

#ifndef SPLITTER_API_H
#define SPLITTER_API_H

#include "core/types_api.h"
#include "core/error_api.h"

/* The <GtSplitter> class defines objects which can split given strings into
   tokens delimited by a given character, allowing for convenient access to
   each token. */
typedef struct GtSplitter GtSplitter;

/* Create a new <GtSplitter> object. */
GtSplitter*   gt_splitter_new(void);

/* Use <splitter> to split <string> of given <length> into tokens delimited by
   <delimiter>. Note that <string> is modified in the splitting process! */
void          gt_splitter_split(GtSplitter *splitter, char *string,
                                GtUword length, char delimiter);
/* Use <splitter> to split <string> of given <length> into tokens delimited by
   <delimiter>. Empty tokens will be ignored. Note that <string> is modified in
   the splitting process! */
void          gt_splitter_split_non_empty(GtSplitter *s, char *string,
                                          GtUword length, char delimiter);
/* Return all tokens split by <splitter> in an array. */
char**        gt_splitter_get_tokens(GtSplitter *splitter);

/* Return token with number <token_num> from <splitter>. */
char*         gt_splitter_get_token(GtSplitter *splitter,
                                    GtUword token_num);
/* Reset the <splitter>. */
void          gt_splitter_reset(GtSplitter *splitter);

/* Return the number of tokens in <splitter>. */
GtUword       gt_splitter_size(GtSplitter *splitter);

/* Delete the <splitter>. */
void          gt_splitter_delete(GtSplitter *splitter);

#endif
