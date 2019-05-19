/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
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

#ifndef MATCH_VISITOR_API_H
#define MATCH_VISITOR_API_H

/* The <GtMatchVisitor> class allows one to distinguish a <GtMatch>
   implementation, e.g. BLAST or OpenMatch, and to call different
   code for each implementation. */
typedef struct GtMatchVisitor GtMatchVisitor;

#include "extended/match_blast_api.h"
#include "extended/match_open_api.h"

/* Visit <match_blast> with <match_visitor>. */
int gt_match_visitor_visit_match_blast(GtMatchVisitor *match_visitor,
                                       GtMatchBlast *match_blast,
                                       GtError *err);

/* Visit <match_open> with <match_visitor>. */
int gt_match_visitor_visit_match_open(GtMatchVisitor *match_visitor,
                                      GtMatchOpen *match_open,
                                      GtError *err);

/* Deletes <match_visitor> freeing all associated space. */
void gt_match_visitor_delete(GtMatchVisitor *match_visitor);

#endif
