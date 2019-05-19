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

#ifndef OUTPUT_FILE_API_H
#define OUTPUT_FILE_API_H

#include "core/file_api.h"
#include "core/option_api.h"

#define GT_FORCE_OPT_CSTR  "force"

/* The <GtOutputFileInfo> class encapsulates output options. */
typedef struct GtOutputFileInfo GtOutputFileInfo;

/* Create a new <GtOutputFileInfo> object. */
GtOutputFileInfo* gt_output_file_info_new(void);
/* Registers the options `-o', `-gzip', `-bzip2' and `-force' in
   <option_parser>. Options chosen during option parsing will be stored in
   <output_file_info> and the output file will be accessible using <*outfp>.
   If no option is given, default <*outfp> will use stdout.
   Caller retains ownership of <*outfp>. */
void              gt_output_file_info_register_options(
                                            GtOutputFileInfo *output_file_info,
                                            GtOptionParser *option_parser,
                                            GtFile **outfp);
/* Deletes <output_file_info> and frees all associated memory. */
void              gt_output_file_info_delete(
                                            GtOutputFileInfo *output_file_info);

#endif
