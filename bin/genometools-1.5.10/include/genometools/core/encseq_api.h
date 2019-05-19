/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCSEQ_API_H
#define ENCSEQ_API_H

#include "core/alphabet_api.h"
#include "core/logger_api.h"
#include "core/timer_api.h"
#include "core/readmode_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"

/* The <GtEncseq> class represents a concatenated collection of sequences from
   one or more input files in a bit-compressed encoding. It is stored in a
   number of <mmap()>able files, depending on which features it is meant to
   support.
   The main compressed sequence information is stored in an __encoded sequence__
   table, with the file suffix '.esq'. This table is the minimum requirement
   for the <GtEncseq> structure and must always be present. In addition, if
   support for multiple sequences is desired, a __sequence separator position__
   table with the '.ssp' suffix is required. If support for sequence
   descriptions is required, two additional tables are needed: a __description__
   table with the suffix '.des' and a __description separator__ table with the
   file suffix '.sds'. Creation and requirement of these tables can be switched
   on and off using API functions as outlined below.
   The <GtEncseq> represents the stored sequences as one concatenated string.
   It allows access to the sequences by providing start positions and lengths
   for each sequence, making it possible to extract encoded substrings into a
   given buffer, as well as accessing single characters both in a random and a
   sequential fashion. */
typedef struct GtEncseq GtEncseq;

/* The <GtEncseqEncoder> class creates objects encapsulating a parameter
   set for conversion from sequence files into encoded sequence files on
   secondary storage. */
typedef struct GtEncseqEncoder GtEncseqEncoder;

/* The <GtEncseqLoader> class creates <GtEncseq> objects by mapping index files
   from secondary storage into memory. */
typedef struct GtEncseqLoader GtEncseqLoader;

/* The <GtEncseqBuilder> class creates <GtEncseq> objects by constructing
   uncompressed, encoded string copies in memory. */
typedef struct GtEncseqBuilder GtEncseqBuilder;

/* The <GtEncseqReader> class represents the current state of a
   sequential scan of a <GtEncseq> region as an iterator. */
typedef struct GtEncseqReader GtEncseqReader;

/* The file suffix used for encoded sequence files. */
#define GT_ENCSEQFILESUFFIX ".esq"
/* The file suffix used for encoded sequence separator position tables. */
#define GT_SSPTABFILESUFFIX ".ssp"
/* The file suffix used for sequence description tables. */
#define GT_DESTABFILESUFFIX ".des"
/* The file suffix used for sequence description separator position tables. */
#define GT_SDSTABFILESUFFIX ".sds"
/* The file suffix used for original input sequence tables. */
#define GT_OISTABFILESUFFIX ".ois"
/* The file suffix used for MD5 fingerprints. */
#define GT_MD5TABFILESUFFIX ".md5"

/* Returns the indexname (as given at loading time) of <encseq> or the string
   "generated" if the GtEncseq was build in memory only. */
const char*       gt_encseq_indexname(const GtEncseq *encseq);
/* Returns the total number of characters in all sequences of <encseq>,
   including separators and wildcards. */
GtUword           gt_encseq_total_length(const GtEncseq *encseq);
/* Returns the total number of sequences contained in <encseq>. */
GtUword           gt_encseq_num_of_sequences(const GtEncseq *encseq);
/* Returns the encoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. */
GtUchar           gt_encseq_get_encoded_char(const GtEncseq *encseq,
                                             GtUword pos,
                                             GtReadmode readmode);
/* Returns the decoded representation of the character at position <pos> of
   <encseq> read in the direction as indicated by <readmode>. */
char              gt_encseq_get_decoded_char(const GtEncseq *encseq,
                                             GtUword pos,
                                             GtReadmode readmode);
/* Returns true iff <pos> is a separator position of <encseq>
   read in the direction as indicated by <readmode>. */
bool              gt_encseq_position_is_separator(const GtEncseq *encseq,
                                                  GtUword pos,
                                                  GtReadmode readmode);
/* Returns true iff <pos> is a wildcard in <encseq>
   read in the direction as indicated by <readmode>. */
bool              gt_encseq_position_is_wildcard(const GtEncseq *encseq,
                                                 GtUword pos,
                                                 GtReadmode readmode);
/* Increases the reference count of <encseq>. */
GtEncseq*         gt_encseq_ref(GtEncseq *encseq);
/* Returns a new <GtEncseqReader> for <encseq>, starting from position
   <startpos>. Also supports reading the sequence from the reverse and
   delivering (reverse) complement characters on DNA alphabets using the
   <readmode> option. Please make sure that the <GT_READMODE_COMPL> and
   <GT_READMODE_REVCOMPL> readmodes are only used on DNA alphabets. */
GtEncseqReader*   gt_encseq_create_reader_with_readmode(const GtEncseq *encseq,
                                                        GtReadmode readmode,
                                                        GtUword startpos);
/* Stores the encoded representation of the substring from 0-based position
   <frompos> to position <topos> of <encseq>. The result is written to the
   location pointed to by <buffer>, which must be large enough to hold the
   result. */
void              gt_encseq_extract_encoded(const GtEncseq *encseq,
                                            GtUchar *buffer,
                                            GtUword frompos,
                                            GtUword topos);
/* Stores the decoded version of the substring from 0-based position <frompos>
   to position <topos> of <encseq>. If the extracted region contains a separator
   character, it will be represented by non-printable SEPARATOR constant.
   The caller is responsible to handle this case. The result of the extraction
   is written to the location pointed to by <buffer>, which must be sufficiently
   large to hold the result. */
void              gt_encseq_extract_decoded(const GtEncseq *encseq,
                                            char *buffer,
                                            GtUword frompos,
                                            GtUword topos);
/* Returns the length of the <seqnum>-th sequence in the <encseq>.
   Requires multiple sequence support enabled in <encseq>. */
GtUword           gt_encseq_seqlength(const GtEncseq *encseq,
                                      GtUword seqnum);
/* Returns the length of the shortest sequence in the <encseq>. */
GtUword           gt_encseq_min_seq_length(const GtEncseq *encseq);
/* Returns the length of the longest sequence in the <encseq>. */
GtUword           gt_encseq_max_seq_length(const GtEncseq *encseq);
/* Returns <true> if <encseq> has multiple sequence support. */
bool              gt_encseq_has_multiseq_support(const GtEncseq *encseq);
/* Returns <true> if <encseq> has description support. */
bool              gt_encseq_has_description_support(const GtEncseq *encseq);
/* Returns <true> if <encseq> has MD5 support. */
bool              gt_encseq_has_md5_support(const GtEncseq *encseq);
/* Returns the start position of the <seqnum>-th sequence in the <encseq>.
   Requires multiple sequence support enabled in <encseq>. */
GtUword           gt_encseq_seqstartpos(const GtEncseq *encseq,
                                        GtUword seqnum);
/* Returns the sequence number from the given <position> for a given
   GtEncseq <encseq>. */
GtUword           gt_encseq_seqnum(const GtEncseq *encseq,
                                   GtUword position);
/* Returns a pointer to the description of the <seqnum>-th sequence in the
   <encseq>. The length of the returned string is written to the
   location pointed at by <desclen>.
   The returned description pointer is not <\0>-terminated!
   Requires description support enabled in <encseq>. */
const char*       gt_encseq_description(const GtEncseq *encseq,
                                        GtUword *desclen,
                                        GtUword seqnum);
/* Returns a <GtStrArray> of the names of the original sequence files
   contained in <encseq>. */
const GtStrArray* gt_encseq_filenames(const GtEncseq *encseq);
/* Returns the number of files contained in <encseq>. */
GtUword           gt_encseq_num_of_files(const GtEncseq *encseq);
/* Returns the effective length (sum of sequence lengths and separators
   between them) of the <filenum>-th file contained in <encseq>. */
GtUint64          gt_encseq_effective_filelength(const GtEncseq *encseq,
                                                 GtUword filenum);
/* Returns the start position of the sequences of the  <filenum>-th file in the
   <encseq>. Requires multiple file support enabled in <encseq>. */
GtUword           gt_encseq_filestartpos(const GtEncseq *encseq,
                                         GtUword filenum);
/* Returns the file number from the given <position> for a given
   GtEncseq <encseq>. */
GtUword           gt_encseq_filenum(const GtEncseq *encseq,
                                    GtUword position);
/* Returns the first sequence number of the sequences in file <filenum> for a
   given GtEncseq <encseq>. */
GtUword           gt_encseq_filenum_first_seqnum(const GtEncseq *encseq,
                                                 GtUword filenum);
/* Returns the <GtAlphabet> associated with <encseq>. */
GtAlphabet*       gt_encseq_alphabet(const GtEncseq *encseq);
/* Extends <encseq>  by virtual reverse complement sequences.
   Returns 0 if mirroring has been successfully enabled, otherwise -1.
   <err> is set accordingly. */
int               gt_encseq_mirror(GtEncseq *encseq, GtError *err);
/* Removes virtual reverse complement sequences added by
   <gt_encseq_mirror()>. */
void              gt_encseq_unmirror(GtEncseq *encseq);
/* Returns <true> if <encseq> contains virtual reverse complement sequences as
   added by <gt_encseq_mirror()>. */
bool              gt_encseq_is_mirrored(const GtEncseq *encseq);
/* Returns the version number of the file representation of <encseq> if it
   exists, or 0 if it was not mapped from a file. */
GtUword           gt_encseq_version(const GtEncseq *encseq);
/* Returns TRUE if <encseq> was created on a 64-bit system. */
bool              gt_encseq_is_64_bit(const GtEncseq *encseq);
/* Deletes <encseq> and frees all associated space. */
void              gt_encseq_delete(GtEncseq *encseq);

/* Reinitializes the given <esr> with the values as described in
   <gt_encseq_create_reader_with_readmode()>. */
void              gt_encseq_reader_reinit_with_readmode(GtEncseqReader *esr,
                                                      const GtEncseq *encseq,
                                                      GtReadmode readmode,
                                                      GtUword startpos);
/* Returns the next encoded character from current position of <esr>, advancing
   the iterator by one position. */
GtUchar           gt_encseq_reader_next_encoded_char(GtEncseqReader *esr);
/* Returns the next decoded character from current position of <esr>, advancing
   the iterator by one position. */
char              gt_encseq_reader_next_decoded_char(GtEncseqReader *esr);
/* Deletes <esr>, freeing all associated space. */
void              gt_encseq_reader_delete(GtEncseqReader *esr);

/* Creates a new <GtEncseqEncoder>. */
GtEncseqEncoder*  gt_encseq_encoder_new(void);
/* Sets <t> to be the timer for <ee>. Default is <NULL> (no progress
reporting). */
void              gt_encseq_encoder_set_timer(GtEncseqEncoder *ee, GtTimer *t);
/* Returns the timer set for <ee>. */
GtTimer*          gt_encseq_encoder_get_timer(const GtEncseqEncoder *ee);
/* Sets the representation of <ee> to <sat> which must be one of 'direct',
   'bytecompress', 'bit', 'uchar', 'ushort' or 'uint32'. Returns 0 on success,
   and a negative value on error (<err> is set accordingly). */
int               gt_encseq_encoder_use_representation(GtEncseqEncoder *ee,
                                                      const char *sat,
                                                      GtError *err);
/* Returns the representation requested for <ee>. */
GtStr*            gt_encseq_encoder_representation(const GtEncseqEncoder *ee);
/* Sets the symbol map file to use in <ee> to <smap> which must a valid
   alphabet description file. Returns 0 on success, and a negative value on
   error (<err> is set accordingly). Default is <NULL> (no alphabet
   transformation). */
int               gt_encseq_encoder_use_symbolmap_file(GtEncseqEncoder *ee,
                                                      const char *smap,
                                                      GtError *err);
/* Returns the symbol map file requested for <ee>. */
const char*       gt_encseq_encoder_symbolmap_file(const GtEncseqEncoder *ee);
/* Sets the logger to use by <ee> during encoding to <l>. Default is <NULL> (no
   logging). */
void              gt_encseq_encoder_set_logger(GtEncseqEncoder *ee,
                                              GtLogger *l);
/* Enables support for retrieving descriptions from the encoded sequence
   encoded by <ee>. That is, the .des and .sds tables are created.
   This is a prerequisite for being able to activate description support in
   <gt_encseq_loader_require_description_support()>. Activated by default. */
void              gt_encseq_encoder_enable_description_support(
                                                           GtEncseqEncoder *ee);
/* Disables support for retrieving descriptions from the encoded sequence
   encoded by <ee>. That is, the .des and .sds tables are not created.
   Encoded sequences created without this support will not be able to be
   loaded via a <GtEncseqLoader> with
   <gt_encseq_loader_require_description_support()> enabled. */
void              gt_encseq_encoder_disable_description_support(
                                                           GtEncseqEncoder *ee);
/* Enables support for random access to multiple sequences in the encoded
   sequence encoded by <ee>. That is, the .ssp table is created.
   This is a prerequisite for being able to activate description support in
   <gt_encseq_loader_require_multiseq_support()>. Activated by default. */
void              gt_encseq_encoder_enable_multiseq_support(
                                                           GtEncseqEncoder *ee);
/* Disables support for random access to multiple sequences in the encoded
   sequence encoded by <ee>. That is, the .ssp table is not created.
   Encoded sequences created without this support will not be able to be
   loaded via a <GtEncseqLoader> with
   <gt_encseq_loader_require_multiseq_support()> enabled. */
void              gt_encseq_encoder_disable_multiseq_support(
                                                           GtEncseqEncoder *ee);
/* Enables support for lossless reproduction of the original sequence,
   regardless of alphabet transformations that may apply. Deactivated by
   default. */
void              gt_encseq_encoder_enable_lossless_support(
                                                           GtEncseqEncoder *ee);
/* Enables support for lossless reproduction of the original sequence,
   regardless of alphabet transformations that may apply. Encoded sequences
   created without this support will not be able to be loaded via a
   <GtEncseqLoader> with <gt_encseq_loader_require_lossless_support()>
   enabled. */
void              gt_encseq_encoder_disable_lossless_support(
                                                           GtEncseqEncoder *ee);
/* Enables support for quick MD5 indexing of the sequences in <ee>. Activated by
   default. */
void              gt_encseq_encoder_enable_md5_support(GtEncseqEncoder *ee);
/* Enables support for quick MD5 indexing of the sequences in <ee>. Encoded
   sequences created without this support will not be able to be loaded via
   a <GtEncseqLoader> with <gt_encseq_loader_require_md5_support()> enabled.*/
void              gt_encseq_encoder_disable_md5_support(GtEncseqEncoder *ee);
/* Enables creation of the .des table containing sequence descriptions.
   Enabled by default. */
void              gt_encseq_encoder_create_des_tab(GtEncseqEncoder *ee);
/* Disables creation of the .des table. */
void              gt_encseq_encoder_do_not_create_des_tab(GtEncseqEncoder *ee);
/* Returns <true> if the creation of the .des table has been requested,
   <false> otherwise. */
bool              gt_encseq_encoder_des_tab_requested(
                                                     const GtEncseqEncoder *ee);
/* Enables creation of the .ssp table containing indexes for multiple sequences.
   Enabled by default. */
void              gt_encseq_encoder_create_ssp_tab(GtEncseqEncoder *ee);
/* Disables creation of the .ssp table. */
void              gt_encseq_encoder_do_not_create_ssp_tab(GtEncseqEncoder *ee);
/* Returns <true> if the creation of the .ssp table has been requested,
   <false> otherwise. */
bool              gt_encseq_encoder_ssp_tab_requested(
                                                     const GtEncseqEncoder *ee);
/* Enables creation of the .sds table containing indexes for sequence
   descriptions. Enabled by default. */
void              gt_encseq_encoder_create_sds_tab(GtEncseqEncoder *ee);
/* Disables creation of the .sds table. */
void              gt_encseq_encoder_do_not_create_sds_tab(GtEncseqEncoder *ee);
/* Returns <true> if the creation of the .sds table has been requested,
   <false> otherwise. */
bool              gt_encseq_encoder_sds_tab_requested(
                                                     const GtEncseqEncoder *ee);
/* Enables creation of the .md5 table containing MD5 sums. Enabled by
   default. */
void              gt_encseq_encoder_create_md5_tab(GtEncseqEncoder *ee);
/* Disables creation of the .md5 table. */
void              gt_encseq_encoder_do_not_create_md5_tab(GtEncseqEncoder *ee);
/* Returns <true> if the creation of the .md5 table has been requested,
   <false> otherwise. */
bool              gt_encseq_encoder_md5_tab_requested(
                                                     const GtEncseqEncoder *ee);
/* Sets the sequence input type for <ee> to DNA. */
void              gt_encseq_encoder_set_input_dna(GtEncseqEncoder *ee);
/* Returns <true> if the input sequence has been defined as being DNA. */
bool              gt_encseq_encoder_is_input_dna(GtEncseqEncoder *ee);
/* Sets the sequence input type for <ee> to protein/amino acids. */
void              gt_encseq_encoder_set_input_protein(GtEncseqEncoder *ee);
/* Returns <true> if the input sequence has been defined as being protein. */
bool              gt_encseq_encoder_is_input_protein(GtEncseqEncoder *ee);
/* Makes <ee> ignore all description suffixes after the first whitespace
   character per description (as defined via isspace(3)). */
void              gt_encseq_encoder_clip_desc(GtEncseqEncoder *ee);
/* Returns <true> if <ee> clips all descriptions after the first whitespace. */
bool              gt_encseq_encoder_are_descs_clipped(GtEncseqEncoder *ee);
/* Encodes the sequence files given in <seqfiles> using the settings in <ee>
   and <indexname> as the prefix for the index tables. Returns 0 on success, or
   a negative value on error (<err> is set accordingly). */
int               gt_encseq_encoder_encode(GtEncseqEncoder *ee,
                                          GtStrArray *seqfiles,
                                          const char *indexname,
                                          GtError *err);
/* Deletes <ee>. */
void              gt_encseq_encoder_delete(GtEncseqEncoder *ee);

/* Creates a new <GtEncseqLoader>. */
GtEncseqLoader*   gt_encseq_loader_new(void);
/* Enables auto-discovery of supported features when loading an encoded
   sequence. That is, if a file with <indexname>.<suffix> exists which
   is named like a table file, it is loaded automatically.
   Use gt_encseq_has_multiseq_support() etc. to query for these capabilities. */
void              gt_encseq_loader_enable_autosupport(GtEncseqLoader *el);
/* Disables auto-discovery of supported features. */
void              gt_encseq_loader_disable_autosupport(GtEncseqLoader *el);
/* Enables support for retrieving descriptions from the encoded sequence
   to be loaded by <el>. That is, the .des and .sds tables must be present.
   For example, these tables are created by having enabled the
   <gt_encseq_encoder_enable_description_support()> option when encoding.
   Activated by default. */
void              gt_encseq_loader_require_description_support(
                                                            GtEncseqLoader *el);
/* Disables support for retrieving descriptions from the encoded sequence
   to be loaded by <el>. That is, the .des and .sds tables need not be present.
   However, disabling this support will result in an error when trying to call
   the method <gt_encseq_description()> on the <GtEncseq>
   object created by <el>. */
void              gt_encseq_loader_drop_description_support(GtEncseqLoader *el);
/* Enables support for random access to multiple sequences in the encoded
   sequence to be loaded by <el>. That is, the .ssp table must be present.
   For example, this table is created by having enabled the
   <gt_encseq_encoder_enable_multiseq_support()> option when encoding.
   Activated by default. */
void              gt_encseq_loader_require_multiseq_support(GtEncseqLoader *el);
/* Disables support for random access to multiple sequences in the encoded
   sequence to be loaded by <el>. That is, the .ssp table needs not be present.
   However, disabling this support will result in an error when trying to call
   the method <gt_encseq_seqlength()> and <gt_encseq_seqstartpos()> on
   the <GtEncseq> object created by <el>. */
void              gt_encseq_loader_drop_multiseq_support(GtEncseqLoader *el);
/* Enables support for lossless reproduction of the original sequence
   in the encoded sequence to be loaded by <el>. That is, the .ois table
   must be present.
   For example, this table is created by having enabled the
   <gt_encseq_encoder_enable_lossless_support()> option when encoding.
   Deactivated by default. */
void              gt_encseq_loader_require_lossless_support(GtEncseqLoader *el);
/* Disables support for lossless reproduction of the original sequence
   in the encoded sequence to be loaded by <el>. That is, the .ois table
   needs not be present.
   However, disabling this support may result in a reduced alphabet
   representation when accessing decoded characters. */
void              gt_encseq_loader_drop_lossless_support(GtEncseqLoader *el);
/* Enables support for quick retrieval of the MD5 sums for the sequences in the
   encoded sequence to be loaded by <el>. That is, the .md5 table must be
   present. For example, this table is created by having enabled the
   <gt_encseq_encoder_enable_md5_support()> option when encoding.
   Activated by default. */
void              gt_encseq_loader_require_md5_support(GtEncseqLoader *el);
/* Disables support for quick retrieval of the MD5 sums for the sequences in the
   encoded sequence to be loaded by <el>. That is, the .md5 table needs not be
   present. */
void              gt_encseq_loader_drop_md5_support(GtEncseqLoader *el);
/* Requires presence of the .des table containing sequence descriptions.
   Enabled by default. */
void              gt_encseq_loader_require_des_tab(GtEncseqLoader *el);
/* Disables requirement of the .des table for loading a <GtEncseq>
   using <el>. */
void              gt_encseq_loader_do_not_require_des_tab(GtEncseqLoader *el);
/* Returns <true> if a .des table must be present for loading to succeed. */
bool              gt_encseq_loader_des_tab_required(const GtEncseqLoader *el);
/* Requires presence of the .ssp table containing indexes for multiple
   sequences. Enabled by default. */
void              gt_encseq_loader_require_ssp_tab(GtEncseqLoader *el);
/* Disables requirement of the .ssp table for loading a <GtEncseq>
   using <el>. */
void              gt_encseq_loader_do_not_require_ssp_tab(GtEncseqLoader *el);
/* Returns <true> if a .ssp table must be present for loading to succeed. */
bool              gt_encseq_loader_ssp_tab_required(const GtEncseqLoader *el);
/* Requires presence of the .sds table containing indexes for sequence
   descriptions. Enabled by default. */
void              gt_encseq_loader_require_sds_tab(GtEncseqLoader *el);
/* Disables requirement of the .sds table for loading a <GtEncseq>
   using <el>. */
void              gt_encseq_loader_do_not_require_sds_tab(GtEncseqLoader *el);
/* Returns <true> if a .sds table must be present for loading to succeed. */
bool              gt_encseq_loader_sds_tab_required(const GtEncseqLoader *el);
/* Sets the logger to use by <ee> during encoding to <l>. Default is <NULL> (no
   logging). */
void              gt_encseq_loader_set_logger(GtEncseqLoader *el, GtLogger *l);
/* Enables loading of a sequence using <el> with mirroring enabled from the
   start. Identical to invoking <gt_encseq_mirror()> directly after loading. */
void              gt_encseq_loader_mirror(GtEncseqLoader *el);
/* Disables loading of a sequence using <el> with mirroring enabled right from
   the start. */
void              gt_encseq_loader_do_not_mirror(GtEncseqLoader *el);
/* Attempts to map the index files as specified by <indexname> using the options
   set in <el> using this interface. Returns a <GtEncseq> instance
   on success, or <NULL> on error. If an error occurred, <err> is set
   accordingly. */
GtEncseq*         gt_encseq_loader_load(GtEncseqLoader *el,
                                        const char *indexname,
                                        GtError *err);
/* Deletes <el>. */
void              gt_encseq_loader_delete(GtEncseqLoader *el);

/* Creates a new <GtEncseqBuilder> using the alphabet <alpha> as a basis for
   on-the-fly encoding of sequences in memory. */
GtEncseqBuilder*  gt_encseq_builder_new(GtAlphabet *alpha);
/* Enables support for retrieving descriptions from the encoded sequence
   to be built by <eb>. Requires additional memory to hold the descriptions and
   a position index.
   Activated by default. */
void              gt_encseq_builder_enable_description_support(
                                                           GtEncseqBuilder *eb);
/* Disables support for retrieving descriptions from the encoded sequence
   to be built by <eb>. Disabling this support will result in an error when
   trying to call the method <gt_encseq_description()> on the
   <GtEncseq> object created by <eb>. */
void              gt_encseq_builder_disable_description_support(
                                                           GtEncseqBuilder *eb);
/* Enables support for random access to multiple sequences in the encoded
   sequence to be built by <eb>. Requires additional memory for an index of
   starting positions. Activated by default. */
void              gt_encseq_builder_enable_multiseq_support(
                                                           GtEncseqBuilder *eb);
/* Disables support for random access to multiple sequences in the encoded
   sequence to be built by <eb>. Disabling this support will result in an
   error when trying to call the method <gt_encseq_seqlength()> or
   <gt_encseq_seqstartpos()> on the <GtEncseq> object created by <eb>. */
void              gt_encseq_builder_disable_multiseq_support(
                                                           GtEncseqBuilder *eb);
/* Enables creation of the .esq table containing the encoded sequence itself.
   Naturally, enabled by default. */
void              gt_encseq_builder_create_esq_tab(GtEncseqBuilder *eb);
/* Disables creation of the .esq table. */
void              gt_encseq_builder_do_not_create_esq_tab(GtEncseqBuilder *eb);
/* Enables creation of the .des table containing sequence descriptions. */
void              gt_encseq_builder_create_des_tab(GtEncseqBuilder *eb);
/* Disables creation of the .des table. */
void              gt_encseq_builder_do_not_create_des_tab(GtEncseqBuilder *eb);
/* Enables creation of the .ssp table containing indexes for multiple sequences.
   */
void              gt_encseq_builder_create_ssp_tab(GtEncseqBuilder *eb);
/* Disables creation of the .ssp table. */
void              gt_encseq_builder_do_not_create_ssp_tab(GtEncseqBuilder *eb);
/* Enables creation of the .sds table containing indexes for sequence
   descriptions. */
void              gt_encseq_builder_create_sds_tab(GtEncseqBuilder *eb);
/* Disables creation of the .sds table. */
void              gt_encseq_builder_do_not_create_sds_tab(GtEncseqBuilder *eb);
/* Adds a sequence given as a C string <str> of length <strlen> to the
   encoded sequence to be built by <eb>. Additionally, a description can be
   given (<desc>). If description support is enabled, this must not be <NULL>.
   A copy will be made during the addition process and the sequence will
   be encoded using the alphabet set at the construction time of <eb>. Thus it
   must only contain symbols compatible with the alphabet. */
void              gt_encseq_builder_add_cstr(GtEncseqBuilder *eb,
                                             const char *str,
                                             GtUword strlen,
                                             const char *desc);
/* Adds a sequence given as a GtStr <str> to the encoded sequence to be built
   by <eb>. Additionally, a description can be given. If description support
   is enabled, <desc> must not be <NULL>.
   A copy will be made during the addition process and the sequence will
   be encoded using the alphabet set at the construction time of <eb>. Thus it
   must only contain symbols compatible with the alphabet. */
void              gt_encseq_builder_add_str(GtEncseqBuilder *eb, GtStr *str,
                                           const char *desc);
/* Adds a sequence given as a pre-encoded string  <str> of length <strlen> to
   the encoded sequence to be built by <eb>. <str> must be encoded using the
   alphabet set at the construction time of <eb>. <str> is not allowed to
   include sequence separators. Does not take ownership of <str>.
   Additionally, a description <desc> can be given. If description support
   is enabled, this must not be <NULL>. */
void              gt_encseq_builder_add_encoded(GtEncseqBuilder *eb,
                                               const GtUchar *str,
                                               GtUword strlen,
                                               const char *desc);
/* Adds a sequence given as a pre-encoded string  <str> of length <strlen> to
   the encoded sequence to be built by <eb>. <str> must be encoded using the
   alphabet set at the construction time of <eb>.
   Always creates a copy of <str>, so it can be used with memory that is to be
   freed immediately after adding.
   Additionally, a description <desc> can be given. If description support
   is enabled, this must not be <NULL>. */
void              gt_encseq_builder_add_encoded_own(GtEncseqBuilder *eb,
                                                   const GtUchar *str,
                                                   GtUword strlen,
                                                   const char *desc);
/* Adds a sequence given as a pre-encoded string  <str> of length <strlen> to
   the encoded sequence to be built by <eb>. <str> must be encoded using the
   alphabet set at the construction time of <eb>. <str> may include
   sequence separators. Does not take ownership of <str>.*/
void              gt_encseq_builder_add_multiple_encoded(GtEncseqBuilder *eb,
                                                         const GtUchar *str,
                                                         GtUword strlen);
/* Sets the logger to use by <ee> during encoding to <l>. Default is <NULL> (no
   logging). */
void              gt_encseq_builder_set_logger(GtEncseqBuilder*, GtLogger *l);

/* Creates a new <GtEncseq> from the sequences added to <eb>.
   Returns a <GtEncseq> instance on success, or <NULL> on error.
   If an error occurred, <err> is set accordingly.
   The state of <eb> is reset to empty after successful creation of a new
   <GtEncseq> (like having called <gt_encseq_builder_reset()>). */
GtEncseq*         gt_encseq_builder_build(GtEncseqBuilder *eb, GtError *err);
/* Clears all added sequences and descriptions, resetting <eb> to a state
   similar to the state immediately after its initial creation.  */
void              gt_encseq_builder_reset(GtEncseqBuilder *eb);
/* Deletes <eb>. */
void              gt_encseq_builder_delete(GtEncseqBuilder *eb);

#endif
