/*
  Copyright (c) 2005-2010 Gordon Gremme <gordon@gremme.org>
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

#ifndef XANSI_API_H
#define XANSI_API_H

#include <stdarg.h>
#include <stdio.h>
#include "core/types_api.h"

/* XANSI module */

/* This module contains wrappers for the functions from the standard library we
   use. These functions always terminate the program if an error occurs.
   That is, one can use this functions without the need to check for errors. */

/* Similar to <atexit(3)>, terminates on error. */
void   gt_xatexit(void (*function)(void));

/* Similar to <fclose(3)>, terminates on error. */
void   gt_xfclose(FILE*);

/* Similar to <fflush(3)>, terminates on error. */
void   gt_xfflush(FILE*);

/* Similar to <fgetc(3)>, terminates on error. */
int    gt_xfgetc(FILE*);

/* Similar to <fgets(3)>, terminates on error. */
char*  gt_xfgets(char *s, int size, FILE *stream);

/* Similar to <fgetpos(3)>, terminates on error. */
void   gt_xfgetpos(FILE*, fpos_t*);

/* Similar to <fopen(3)>, terminates on error. */
FILE*  gt_xfopen(const char *path, const char *mode);

/* Similar to <fputc(3)>, terminates on error. */
void   gt_xfputc(int, FILE*);

/* Similar to <fputs(3)>, terminates on error. */
void   gt_xfputs(const char*, FILE*);

/* Similar to <fread(3)>, terminates on error. */
size_t gt_xfread(void *ptr, size_t size, size_t nmemb, FILE *fp);

/* Shortcut to <gt_xfread()> which reads a single element of data (of size
   <sizeof (*ptr)>) from <fp> and stores the result in <ptr>. */
#define gt_xfread_one(ptr, fp)                  \
        gt_xfread(ptr, sizeof (*ptr), (size_t) 1, fp)

/* Similar to <fseek(3)>, terminates on error. */
void   gt_xfseek(FILE*, GtWord offset, int whence);

/* Similar to <fsetpos(3)>, terminates on error. */
void   gt_xfsetpos(FILE*, const fpos_t*);

/* Similar to <fwrite(3)>, terminates on error. */
void   gt_xfwrite(const void *ptr, size_t size, size_t nmemb, FILE *fp);

/* Shortcut to <gt_xfwrite()> which writes a single element of data (of size
   <sizeof (*ptr)>) from <ptr> to <fp>. */
#define gt_xfwrite_one(ptr, fp)                  \
        gt_xfwrite(ptr, sizeof (*ptr), (size_t) 1, fp)

/* Similar to <putchar(3)>, terminates on error. */
void   gt_xputchar(int);

/* Similar to <puts(3)>, terminates on error. */
void   gt_xputs(const char*);

/* Similar to <remove(3)>, terminates on error. */
void   gt_xremove(const char*);

/* Similar to <ungetc(3)>, terminates on error. */
void   gt_xungetc(int, FILE*);

/* Similar to <vfprintf(3)>, terminates on error. */
void   gt_xvfprintf(FILE *stream, const char *format, va_list ap);

/* Similar to <vsnprintf(3)>, terminates on error. */
int    gt_xvsnprintf(char *str, size_t size, const char *format, va_list ap);

#endif
