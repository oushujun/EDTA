/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ALPHABET_API_H
#define ALPHABET_API_H

#include <limits.h>
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "core/types_api.h"
#include "core/error_api.h"

/* The following type is for storing alphabets.*/
typedef struct GtAlphabet GtAlphabet;

/* Return a <GtAlphabet> object which represents a DNA alphabet. */
GtAlphabet*    gt_alphabet_new_dna(void);
/* Return a <GtAlphabet> object which represents a protein alphabet. */
GtAlphabet*    gt_alphabet_new_protein(void);
/* Return an empty <GtAlphabet> object. */
GtAlphabet*    gt_alphabet_new_empty(void);
/* Return a <GtAlphabet> object, as read from an .al1 file specified by
   <filename> (i.e. no al1 suffix necessary). */
GtAlphabet*    gt_alphabet_new_from_file(const char *filename, GtError *err);
/* Return a <GtAlphabet> object, as read from a file specified by
   <filename>. */
GtAlphabet*    gt_alphabet_new_from_file_no_suffix(const char *filename,
                                                   GtError *err);
/* Return a <GtAlphabet> object, as read from a string of length <len>
   specified by <alphadef>. */
GtAlphabet*    gt_alphabet_new_from_string(const char *alphadef,
                                           GtUword len, GtError *err);
/* Returns a new <GtAlphabet> object by scanning the sequence files in
   <filenametab> to determine whether they are DNA or protein sequences,
   and the appropriate alphabet will be used (see gt_alphabet_guess()).
   Returns NULL on error, see <err> for details. */
GtAlphabet*    gt_alphabet_new_from_sequence(const GtStrArray *filenametab,
                                             GtError *err);
/* Try to guess which type the given <sequence> with <length> has (DNA or
   protein) and return an according <GtAlphabet*> object. */
GtAlphabet*    gt_alphabet_guess(const char *sequence, GtUword seqlen);
/* Return a clone of <alphabet>. */
GtAlphabet*    gt_alphabet_clone(const GtAlphabet *alphabet);
/* Returns TRUE if <a> and <b> are equal (i.e. have the same symbol
   mapping), FALSE otherwise. */
bool           gt_alphabet_equals(const GtAlphabet *a, const GtAlphabet *b);
/* Increase the reference count for <alphabet> and return it. */
GtAlphabet*    gt_alphabet_ref(GtAlphabet *alphabet);
/* Add the mapping of all given <characters> to the given <alphabet>. The first
   character is the result of subsequent <gt_alphabet_decode()> calls. */
void           gt_alphabet_add_mapping(GtAlphabet *alphabet,
                                       const char *characters);
/* Add <wildcard> to the <alphabet>. */
void           gt_alphabet_add_wildcard(GtAlphabet *alphabet, char wildcard);
/* Returns the array of symbols from <alphabet> such that the index of the
   character equals its encoding. */
const GtUchar* gt_alphabet_symbolmap(const GtAlphabet *alphabet);
/* Returns number of characters in <alphabet> (excluding wildcards). */
unsigned int   gt_alphabet_num_of_chars(const GtAlphabet *alphabet);
/* Returns number of characters in <alphabet> (including wildcards). */
unsigned int   gt_alphabet_size(const GtAlphabet *alphabet);
/* Returns an array of the characters in <alphabet>. */
const GtUchar* gt_alphabet_characters(const GtAlphabet *alphabet);
/* Returns the character used in <alphabet> to represent wildcards in output. */
GtUchar        gt_alphabet_wildcard_show(const GtAlphabet *alphabet);
/* Returns the required number of bits required to represent a symbol
   in <alphabet>. */
unsigned int   gt_alphabet_bits_per_symbol(const GtAlphabet *alphabet);
/* Writes a representation of <alphabet> to the file pointer <fpout>. */
void           gt_alphabet_output(const GtAlphabet *alphabet, FILE *fpout);
/* Writes a representation of <alphabet> to the .al1 output file as specified
   by <indexname> (i.e. without the .al1 suffix). */
int            gt_alphabet_to_file(const GtAlphabet *alphabet,
                                   const char *indexname,
                                   GtError *err);
/* Writes a representation of <alphabet> to the <GtStr> as specified
   by <dest>. */
void           gt_alphabet_to_str(const GtAlphabet *alphabet, GtStr *dest);
/* Returns the printable character specified in <alphabet> for <currentchar>. */
GtUchar        gt_alphabet_pretty_symbol(const GtAlphabet *alphabet,
                                         unsigned int currentchar);
/* Prints the printable character specified in <alphabet> for <currentchar> on
   <fpout>. */
void           gt_alphabet_echo_pretty_symbol(const GtAlphabet *alphabet,
                                              FILE *fpout,
                                              GtUchar currentchar);
/* The following method checks if the given <alphabet> is the protein
   alphabet with the aminoacids A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S,
   T, V, W, Y written in lower or upper case and returns <true>, if this is the
   case (<false> otherwise). */
bool           gt_alphabet_is_protein(const GtAlphabet *alphabet);
/* The following method checks if the given alphabet is the DNA alphabet with
   the bases A, C, G, T written in lower or upper case and returns <true>, if
   this is the case (<false> otherwise). */
bool           gt_alphabet_is_dna(const GtAlphabet *alphabet);
/* Returns true if the character <c> is defined in <alphabet>. */
bool           gt_alphabet_valid_input(const GtAlphabet *alphabet, char c);
/* Encode character <c> with given <alphabet>.
   Ensure that <c> is encodable with the given <alphabet>! */
GtUchar        gt_alphabet_encode(const GtAlphabet *alphabet, char c);
/* Decode character <c> with given <alphabet>. */
char           gt_alphabet_decode(const GtAlphabet *alphabet, GtUchar c);
/* Encode sequence <in> of given <length> with <alphabet> and store the result
   in <out>. <in> has to be encodable with the given <alphabet>! */
void           gt_alphabet_encode_seq(const GtAlphabet *alphabet, GtUchar *out,
                                      const char *in, GtUword length);
/* Suppose the string <src> of length <len> was transformed according to the
   <alphabet>. The following method shows each character in <src> as the
   printable character specified in the transformation. The output is written
   to the given file pointer <fpout>. */
void           gt_alphabet_decode_seq_to_fp(const GtAlphabet *alphabet,
                                            FILE *fpout,
                                            const GtUchar *src,
                                            GtUword len);
/* Suppose the string <src> of length <len> was transformed according to the
   <alphabet>. The following method shows each character in <src> as the
   printable character specified in the transformation. The output is written
   to the given string <dest> and terminated with '\0' at dest[len]. <dest>
   therefore has to be of at least <len> + 1 length. */
void           gt_alphabet_decode_seq_to_cstr(const GtAlphabet *alphabet,
                                              char *dest,
                                              const GtUchar *src,
                                              GtUword len);
/* Analog to <gt_alphabet_decode_seq_to_fp()> writing the output to
   a new <GtStr>. */
GtStr*         gt_alphabet_decode_seq_to_str(const GtAlphabet *alphabet,
                                             const GtUchar *src,
                                             GtUword len);
/* Decrease the reference count for <alphabet> or delete it, if this was the
   last reference. */
void           gt_alphabet_delete(GtAlphabet *alphabet);

#endif
