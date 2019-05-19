/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef TOOLBOX_API_H
#define TOOLBOX_API_H

#include "core/error_api.h"
#include "core/tool_api.h"

/* The <GtToolbox> class groups several tools into one and can be used to
   structure __GenomeTools__ into sensible sets of subtools. */
typedef struct GtToolbox GtToolbox;

/* Return a new empty <GtToolbox>. */
GtToolbox* gt_toolbox_new(void);
/* Add <tool> with name <toolname> to <toolbox>. Takes ownership of <tool>. */
void       gt_toolbox_add_tool(GtToolbox *toolbox, const char *toolname,
                               GtTool *tool);
/* Add (hidden) <tool> with name <toolname> to <toolbox>. Hidden tools are not
   shown in the output of <gt_toolbox_show()>. Takes ownership of <tool>. */
void       gt_toolbox_add_hidden_tool(GtToolbox *toolbox, const char *toolname,
                                      GtTool *tool);
/* Get <GtTool> with name <toolname> from <toolbox>. Returns NULL if tool does
   not exist in <toolbox>. */
GtTool*    gt_toolbox_get_tool(GtToolbox *toolbox, const char *toolname);
/* Show all tools in <toolbox> except the hidden ones. Intended to be used
   as an argument to <gt_option_parser_set_comment_func()>. */
int        gt_toolbox_show(const char *progname, void *toolbox, GtError*);
/* Deletes <toolbox>. */
void       gt_toolbox_delete(GtToolbox *toolbox);

#endif
