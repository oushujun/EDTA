/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef FILEUTILS_API_H
#define FILEUTILS_API_H

#include <stdbool.h>
#include <stdio.h>
#include <sys/types.h>
#include "core/str_api.h"
#include "core/str_array_api.h"

/* Fileutils module */

/* Returns the suffix of <path>, if there is any. Returns "" otherwise.
   The suffix is the part after and including the last '.' but after the last
   '/' (or '\' on Windows). Except if <path> ends with ".gz" or ".bz2", then the
   suffix is the part after and including the second last '.'. */
const char* gt_file_suffix(const char *path);

/* Returns true if the file with the given <path> exists, false otherwise. */
bool        gt_file_exists(const char *path);

/* Returns true if the file with the name composed of the concatenation of
   <path> and <suffix> exists, false otherwise. */
bool        gt_file_exists_with_suffix(const char *path, const char *suffix);

/* Returns true if the file with path <a> has a later modification time than the
   file with path <b>, false otherwise. */
bool        gt_file_is_newer(const char *a, const char *b);

/* Returns the number of lines in a file. */
GtUword     gt_file_number_of_lines(const char*);

/* Set <path> to the dirname of <file>, if it has one, to "" otherwise. */
void        gt_file_dirname(GtStr *path, const char *file);

/* Find <file> in $PATH, if it has no dirname; set <path> to dirname otherwise.
   Sets <path> to the empty string if <file> could not be found in $PATH. */
int         gt_file_find_in_path(GtStr *path, const char *file, GtError*);

/* Find  <file> in the ':'-separated directory list (on Windows ';'-separated)
   specified in environment variable $<env>, if it has no dirname; set <path> to
   dirname otherwise. Sets <path> to the empty string if <file> could not be
   found in $<env>. */
int         gt_file_find_in_env(GtStr *path, const char *file, const char *env,
                                GtError*);

/* Return the (estimated) size of <file>. If <file> is uncompressed, the exact
   size is returned. If <file> is compressed, an estimation which assumes that
   <file> contains a DNA sequence is returned. */
off_t       gt_file_estimate_size(const char *file);

/* Return the (estimated) total size of all files given in <filenames>.
   Uses <gt_file_estimate_size()>. */
off_t       gt_files_estimate_total_size(const GtStrArray *filenames);

/* Guesse if the sequences contained in the files given in <filenames> are
   protein sequences. Returns 1 if the guess is that the files contain protein
   sequences. Returns 0 if the guess is that the files contain DNA sequences.
   Returns -1 if an error occurs while reading the files (<err> is set
   accordingly). */
int         gt_files_guess_if_protein_sequences(const GtStrArray *filenames,
                                                GtError *err);

/* Find regular executable <file> in $PATH, if it has no dirname; set <path> to
   dirname otherwise. Sets <path> to the empty string if regular executable
   <file> could not be found in $PATH. */
int         gt_file_find_exec_in_path(GtStr *path, const char *file,
                                      GtError *err);

/* Return the size of <file> in bytes. <file> must exist. */
off_t       gt_file_size(const char *file);

/* Returns the size of the file whose name is composed of the
  concatenation of <path> and <suffix>. <file> must exist. */
off_t       gt_file_size_with_suffix(const char *path, const char *suffix);

#endif
