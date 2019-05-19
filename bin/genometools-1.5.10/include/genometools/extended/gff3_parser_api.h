/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef GFF3_PARSER_API_H
#define GFF3_PARSER_API_H

#include "core/cstr_table_api.h"
#include "core/file_api.h"
#include "core/queue_api.h"
#include "core/range_api.h"
#include "core/strand_api.h"
#include "core/types_api.h"
#include "extended/type_checker_api.h"
#include "extended/xrf_checker_api.h"

/* A <GtGFF3Parser> can be used to parse GFF3 files and convert them into
   <GtGenomeNode> objects. If the GFF3 files do not contain the encouraged
   sequence-region meta directives, the GFF3 parser introduces the corresponding
   region nodes automatically. This is a low-level class and it is usually not
   used directly. Normally, a <GtGFF3InStream> is used to parse GFF3 files. */
typedef struct GtGFF3Parser GtGFF3Parser;

/* Return a new <GtGFF3Parser> object with optional <type_checker>. If a
   <type_checker> was given, the <GtGFF3Parser> stores a new reference to it
   internally and uses the <type_checker> to check types during parsing. */
GtGFF3Parser* gt_gff3_parser_new(GtTypeChecker *type_checker);
/* Enable ID attribute checking in <gff3_parser>. Thereby, the memory
   consumption of the <gff3_parser> becomes proportional to the input file
   size(s). */
void          gt_gff3_parser_check_id_attributes(GtGFF3Parser *gff3_parser);
/* Enable sequence region boundary checking in <gff3_parser>. That is,
   encountering features outside the sequence region boundaries will result in
   an error. */
void          gt_gff3_parser_check_region_boundaries(GtGFF3Parser *gff3_parser);
/* Disable sequence region boundary checking in <gff3_parser>. That is,
   features outside the sequence region boundaries will be permitted. */
void          gt_gff3_parser_do_not_check_region_boundaries(GtGFF3Parser
                                                                  *gff3_parser);
/* Transform all features parsed by <gff3_parser> by the given <offset>. */
void          gt_gff3_parser_set_offset(GtGFF3Parser *gff3_parser,
                                        GtWord offset);
/* Set <type_checker> used by <gff3_parser>. */
void          gt_gff3_parser_set_type_checker(GtGFF3Parser *gff3_parser,
                                              GtTypeChecker *type_checker);
/* Set <xrf_checker> used by <gff3_parser>. */
void          gt_gff3_parser_set_xrf_checker(GtGFF3Parser *gff3_parser,
                                             GtXRFChecker *xrf_checker);
/* Enable the tidy mode in <gff3_parser>. In tidy mode the <gff3_parser> parser
   tries to tidy up features which would normally lead to a parse error. */
void          gt_gff3_parser_enable_tidy_mode(GtGFF3Parser *gff3_parser);
/* Use <gff3_parser> to parse genome nodes from file pointer <fpin>.
   <status_code> is set to 0 if at least one genome node was created (and stored
   in <genome_nodes>) and to <EOF> if no further genome nodes could be parsed
   from <fpin>. Every encountered (genome feature) type is recorded in the
   C string table <used_types>. The parser uses the given <filenamestr> to
   store the file name of <fpin> in the created genome nodes or to give the
   correct filename in error messages, if necessary.
   <line_number> is increased accordingly during parsing and has to be set to 0
   before parsing a new <fpin>.
   If an error occurs during parsing this method returns -1 and sets <err>
   accordingly. */
int           gt_gff3_parser_parse_genome_nodes(GtGFF3Parser *gff3_parser,
                                                int *status_code,
                                                GtQueue *genome_nodes,
                                                GtCstrTable *used_types,
                                                GtStr *filenamestr,
                                                GtUint64 *line_number,
                                                GtFile *fpin,
                                                GtError *err);
/* Reset the <gff3_parser> (necessary if the input file is switched). */
void          gt_gff3_parser_reset(GtGFF3Parser *gff3_parser);
/* Delete the <gff3_parser>. */
void          gt_gff3_parser_delete(GtGFF3Parser *gff3_parser);

#endif
