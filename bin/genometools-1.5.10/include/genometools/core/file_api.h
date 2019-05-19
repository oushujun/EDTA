/*
  Copyright (c) 2005-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef FILE_API_H
#define FILE_API_H

#include <stdio.h>
#include "core/error_api.h"

/* This class defines (generic) files in __GenomeTools__. A generic file is is a
   file which either uncompressed or compressed (with gzip or bzip2).
   A <NULL>-pointer as generic file implies <stdout>. */
typedef struct GtFile GtFile;

/* Return a new <GtFile> object for the given <path> and open the underlying
   file handle with given <mode>. Returns <NULL> and sets <err> accordingly, if
   the file <path> could not be opened. The compression mode is determined by
   the ending of <path> (gzip compression if it ends with '.gz', bzip2
   compression if it ends with '.bz2', and uncompressed otherwise). */
GtFile* gt_file_new(const char *path, const char *mode, GtError *err);

/* Increments the reference count of <file>. */
GtFile* gt_file_ref(GtFile *file);

/* Create a new <GtFile> object from a normal file pointer <fp>. */
GtFile* gt_file_new_from_fileptr(FILE *fp);

/* <printf(3)> for generic <file>. */
void    gt_file_xprintf(GtFile *file, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));

/* Write <\0>-terminated C string <cstr> to <file>. Similar to <fputs(3)>, but
   terminates on error. */
void    gt_file_xfputs(const char *cstr, GtFile *file);

/* Write single character <c> to <file>. Similar to <fputc(3)>, but terminates
   on error. */
void    gt_file_xfputc(int c, GtFile *file);

/* Return next character from <file> or <EOF>, if end-of-file is reached. */
int     gt_file_xfgetc(GtFile *file);

/* Read up to <nbytes> from generic <file> and store result in <buf>, returns
   bytes read. */
int     gt_file_xread(GtFile *file, void *buf, size_t nbytes);

/* Write <nbytes> from <buf> to given generic <file>. */
void    gt_file_xwrite(GtFile *file, void *buf, size_t nbytes);

/* Rewind the generic <file>. */
void    gt_file_xrewind(GtFile *file);

/* Close the underlying file handle and destroy the <file> object. */
void    gt_file_delete(GtFile *file);

/* Destroy the file handle object, but do not close the underlying handle. */
void    gt_file_delete_without_handle(GtFile*);

#endif
