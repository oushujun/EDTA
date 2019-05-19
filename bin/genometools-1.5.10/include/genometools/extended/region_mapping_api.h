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

#ifndef REGION_MAPPING_API_H
#define REGION_MAPPING_API_H

#include "core/encseq_api.h"
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"

/* A <GtRegionMapping> objects maps sequence-regions to the corresponding
   entries of sequence files. */
typedef struct GtRegionMapping GtRegionMapping;

/* Return a new <GtRegionMapping> object for the mapping file with the given
   <mapping_filename>. In the case of an error, <NULL> is returned and <err> is
   set accordingly. */
GtRegionMapping* gt_region_mapping_new_mapping(GtStr *mapping_filename,
                                               GtError *err);

/* Return a new <GtRegionMapping> object for the sequence files given in
   <sequence_filenames>. If <matchdesc> is <true>, the sequence descriptions
   from the input files are matched for the desired sequence IDs (in GFF3).

   If <usedesc> is <true>, the sequence descriptions are used to map the
   sequence IDs (in GFF3) to actual sequence entries. If a description contains
   a sequence range (e.g., III:1000001..2000000), the first part is used as
   sequence ID ('III') and the first range position as offset ('1000001').

   <matchdesc> and <usedesc> cannot be <true> at the same time. */
GtRegionMapping* gt_region_mapping_new_seqfiles(GtStrArray *sequence_filenames,
                                                bool matchdesc, bool usedesc);

void             gt_region_mapping_match_start(GtRegionMapping *region_mapping);

/* Like <gt_region_mapping_new_seqfiles()>, but using <encseq> as a sequence
   source. */
GtRegionMapping* gt_region_mapping_new_encseq(GtEncseq *encseq, bool matchdesc,
                                              bool usedesc);

/* Return a new <GtRegionMapping> object which maps to the given sequence
   <rawseq> with the corresponding <length> and <offset>. */
GtRegionMapping* gt_region_mapping_new_rawseq(const char *rawseq,
                                              GtUword length, GtUword offset);

/* Increase the reference count for <region_mapping> and return it. */
GtRegionMapping* gt_region_mapping_ref(GtRegionMapping *region_mapping);

/* Use <region_mapping> to extract the sequence from <start> to <end> of the
   given sequence ID <seqid> into a buffer written to <seq> (the caller is
   responsible to free it).
   In the case of an error, -1 is returned and <err> is set accordingly. */
int              gt_region_mapping_get_sequence(GtRegionMapping *region_mapping,
                                                char **seq, GtStr *seqid,
                                                GtUword start,
                                                GtUword end,
                                                GtError *err);

/* Use <region_mapping> to retrieve the sequence length of the given
   sequence ID <seqid> and store the result in <length>.
   In the case of an error, -1 is returned and <err> is set accordingly. */
int              gt_region_mapping_get_sequence_length(GtRegionMapping
                                                       *region_mapping,
                                                       GtUword *length,
                                                       GtStr *seqid,
                                                       GtError *err);

/* Use <region_mapping> to get the description of the MD5 sequence ID <seqid>.
   The description is appended to <desc>.
   In the case of an error, -1 is returned and <err> is set accordingly. */
int              gt_region_mapping_get_description(GtRegionMapping
                                                   *region_mapping,
                                                   GtStr *desc,
                                                   GtStr *seqid,
                                                   GtError *err);

/* Use <region_mapping> to return the MD5 fingerprint of the sequence with the
   sequence ID <seqid> and its corresponding <range>. The offset of the sequence
   is stored in <offset>.
   In the case of an error, <NULL> is returned and <err> is set accordingly. */
const char*      gt_region_mapping_get_md5_fingerprint(GtRegionMapping
                                                       *region_mapping,
                                                       GtStr *seqid,
                                                       const GtRange *range,
                                                       GtUword *offset,
                                                       GtError *err);

/* Delete <region_mapping>. */
void             gt_region_mapping_delete(GtRegionMapping *region_mapping);

#endif
