/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef MD5_ENCODER_API_H
#define MD5_ENCODER_API_H

/* The <GtMD5Encoder> class implements a stateful encoder for MD5 hashes
   for strings build by iterative addition of blocks. */
typedef struct GtMD5Encoder GtMD5Encoder;

/* Returns a new <GtMD5Encoder> object. */
GtMD5Encoder* gt_md5_encoder_new(void);
/* Processes <message> of length <len> (max. block length 64 bytes) to be
   incorporated in the hash currently represented by <enc>. */
void          gt_md5_encoder_add_block(GtMD5Encoder *enc, const char *message,
                                       GtUword len);
/* Finishes <enc> to produce the final MD5 value, written to the 16-byte array
   <output>. If <outstr> is not NULL, a \0-terminated string representation of
   the hash will be written to the 32-byte string buffer it points to. */
void          gt_md5_encoder_finish(GtMD5Encoder *enc, unsigned char *output,
                                    char *outstr);
/* Discards the current state of <enc> and resets it to represent the MD5 hash
   of the empty string. */
void          gt_md5_encoder_reset(GtMD5Encoder *enc);
/* Deletes <enc> and frees all associated space. */
void          gt_md5_encoder_delete(GtMD5Encoder *enc);
#endif
