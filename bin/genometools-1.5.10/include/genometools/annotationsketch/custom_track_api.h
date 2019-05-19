/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CUSTOM_TRACK_API_H
#define CUSTOM_TRACK_API_H

/* The <GtCustomTrack> interface allows the <GtCanvas> to call user-defined
   drawing functions on a <GtGraphics> object. Please refer to the specific
   implementations' documentation for more information on a particular
   custom track. */
typedef struct GtCustomTrack GtCustomTrack;

/* Increase the reference count for <ctrack>. */
GtCustomTrack* gt_custom_track_ref(GtCustomTrack *ctrack);
/* Delete the given <ctrack>. */
void           gt_custom_track_delete(GtCustomTrack *ctrack);

#endif
