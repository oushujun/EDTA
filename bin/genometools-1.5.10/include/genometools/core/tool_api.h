/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef TOOL_API_H
#define TOOL_API_H

#include "core/error_api.h"
#include "core/option_api.h"

/* The <GtTool> class encapsulates a single __GenomeTools__ tool. Can also be
   used in external applications based on <libgenometools>. */
typedef struct GtTool GtTool;

/* Callback function. Must return memory to be used as a storage space
   for tool arguments. */
typedef void*           (*GtToolArgumentsNew)(void);
/* Callback function. Must free up all memory reserved by the
   <GtToolArgumentsNew> function in <tool_arguments>. */
typedef void            (*GtToolArgumentsDelete)(void *tool_arguments);
/* Callback function. Must return a new <GtOptionParser> filling
   <tool_arguments> with content. */
typedef GtOptionParser* (*GtToolOptionParserNew)(void *tool_arguments);
/* Callback function. Checks the validity of <tool_arguments> when <rest_argc>
   additional parameters are given to the tool command line.
   Must return zero if checks are successful, and a negative value otherwise.
   In that case <err> should be set accordingly. */
typedef int             (*GtToolArgumentsCheck)(int rest_argc,
                                                void *tool_arguments,
                                                GtError *err);
/* Callback function. Acts as a main entry point for the tool logic.
   Parameters <argc> and <argv> are similar to a regular main() function.
   The <parsed_args> parameter gives the number of parameters already parsed
   by the option parser before additional parameters start. Use <tool_arguments>
   to access options set by the option parser and write errors to <err>.
   This function should return the error status of the tool (i.e. 0 for
   success). If the return value is not equal to 0, errors written to <err>
   will be printed on stderr. */
typedef int             (*GtToolRunner)(int argc, const char **argv,
                                        int parsed_args, void *tool_arguments,
                                        GtError*);
/* Returns a new self-contained <GtTool>. */
typedef GtTool*         (*GtToolConstructor)(void);

/* Create a new tool object, with
   a tool argument constructor <gt_tool_arguments_new> (optional),
   a tool argument destructor <gt_tool_arguments_delete> (optional).
   a tool option parser constructor <gt_tool_option_parser_new> (required),
   a tool argument checker <gt_tool_arguments_check> (optional),
   a tool runner <gt_tool_runner> (required), and
   <tool_arguments_new> and <tool_arguments_delete> imply each other.
   Returns a new GtTool object. */
GtTool* gt_tool_new(GtToolArgumentsNew tool_arguments_new,
                    GtToolArgumentsDelete tool_arguments_delete,
                    GtToolOptionParserNew tool_option_parser_new,
                    GtToolArgumentsCheck tool_arguments_check,
                    GtToolRunner tool_runner);

/* Run the given <tool> as follows:
   1. Create a tool arguments object, if necessary.
   2. Create a new option parser and pass the tool arguments along.
   3. Parse the options (<argc> and <argv>) with the created option parser.
      Return upon error, continue otherwise.
   4. Check the tool arguments, if necessary.
      Return upon error, continue otherwise.
   5. Run the actual tool with the given arguments (the tool arguments object
      is passed along).
   6. Delete the tool arguments object, if one was created.
      Returns -1 and sets <err> on error, returns 0 otherwise. */
int     gt_tool_run(GtTool*, int argc, const char **argv, GtError *err);

/* Delete the given <tool>. */
void    gt_tool_delete(GtTool*);

#endif
