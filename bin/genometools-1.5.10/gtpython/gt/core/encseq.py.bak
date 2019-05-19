#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

from gt.dlload import gtlib
from gt.core.alphabet import Alphabet
from gt.core.error import Error, gterror
from gt.core.readmode import *
from gt.core.str_array import StrArray
from gt.core.gtstr import Str
import os
from ctypes import Structure, c_ulong, c_uint, c_int, c_char_p, c_void_p, \
    POINTER, byref, c_uint32, c_uint64, string_at, \
    create_string_buffer, c_ubyte, c_char


def int2bool(val):
    if val == 0:
        return False
    else:
        return True


class Seqinfo(Structure):
    _fields_ = [("startpos", c_ulong), ("length", c_ulong)]


class EncseqEncoder:

    def __init__(self):
        self.ee = gtlib.gt_encseq_encoder_new()
        self._as_parameter_ = self.ee

    def __del__(self):
        try:
            gtlib.gt_encseq_encoder_delete(self.ee)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, EncseqEncoder):
            raise TypeError("argument must be an EncseqEncoder object")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def set_representation_type(self, sat):
        gtlib.gt_encseq_encoder_set_access_type(
            self.ee, str(sat.encode("utf-8")))

    def set_symbolmap_file(self, fn):
        gtlib.gt_encseq_encoder_use_symbolmap_file(
            self.ee, str(fn.encode("utf-8")))

    def enable_description_support(self):
        gtlib.gt_encseq_encoder_enable_description_support(self.ee)

    def disable_description_support(self):
        gtlib.gt_encseq_encoder_disable_description_support(self.ee)

    def enable_multiseq_support(self):
        gtlib.gt_encseq_encoder_enable_multiseq_support(self.ee)

    def disable_multiseq_support(self):
        gtlib.gt_encseq_encoder_disable_multiseq_support(self.ee)

    def create_des_tab(self):
        gtlib.gt_encseq_encoder_create_des_tab(self.ee)

    def do_not_create_des_tab(self):
        gtlib.gt_encseq_encoder_do_not_create_des_tab(self.ee)

    def create_ssp_tab(self):
        gtlib.gt_encseq_encoder_create_ssp_tab(self.ee)

    def do_not_create_ssp_tab(self):
        gtlib.gt_encseq_encoder_do_not_create_ssp_tab(self.ee)

    def create_sds_tab(self):
        gtlib.gt_encseq_encoder_create_sds_tab(self.ee)

    def do_not_create_sds_tab(self):
        gtlib.gt_encseq_encoder_do_not_create_sds_tab(self.ee)

    def set_input_protein(self):
        gtlib.gt_encseq_encoder_set_input_protein(self.ee)

    def set_input_dna(self):
        gtlib.gt_encseq_encoder_set_input_dna(self.ee)

    def encode(self, seqfiles, indexname):
        if not isinstance(seqfiles, list):
            raise TypeError("argument must be a list of strings")
        if len(seqfiles) == 0:
            raise Error("list of input sequence files must be non-empty")
        sa = StrArray()
        for f in seqfiles:
            if not os.path.exists(str(f)):
                raise IOError("file not found: %s" % str(f))
            sa.add(str(f))
        err = Error()
        esptr = gtlib.gt_encseq_encoder_encode(self.ee, sa.strarr,
                                               str(indexname).encode("UTF-8"),
                                               err.error)
        if esptr != 0:
            gterror(err)

    def register(cls, gtlib):
        gtlib.gt_encseq_encoder_new.restype = c_void_p
        gtlib.gt_encseq_encoder_encode.argtypes = [c_void_p, c_void_p,
                                                   c_char_p, c_void_p]
        gtlib.gt_encseq_encoder_encode.restype = c_int
        gtlib.gt_encseq_encoder_delete.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_use_representation.argtypes = [
            c_void_p, c_char_p]
        gtlib.gt_encseq_encoder_use_symbolmap_file.argtypes = [
            c_void_p, c_char_p]
        gtlib.gt_encseq_encoder_enable_description_support.argtypes = [
            c_void_p]
        gtlib.gt_encseq_encoder_disable_description_support.argtypes = [
            c_void_p]
        gtlib.gt_encseq_encoder_enable_multiseq_support.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_disable_multiseq_support.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_create_des_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_do_not_create_des_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_create_ssp_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_do_not_create_ssp_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_create_sds_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_do_not_create_sds_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_set_input_protein.argtypes = [c_void_p]
        gtlib.gt_encseq_encoder_set_input_dna.argtypes = [c_void_p]

    register = classmethod(register)


class EncseqLoader:

    def __init__(self):
        self.el = gtlib.gt_encseq_loader_new()
        self._as_parameter_ = self.el
        self.destab = True
        self.ssptab = True
        self.sdstab = True

    def __del__(self):
        try:
            pass  # FIXME gtlib.gt_encseq_loader_delete(self.el)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, EncseqLoader):
            raise TypeError("argument must be an EncseqLoader object")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def require_description_support(self):
        self.destab = True
        self.sdstab = True
        gtlib.gt_encseq_loader_require_description_support(self.el)

    def drop_description_support(self):
        self.destab = False
        self.sdstab = False
        gtlib.gt_encseq_loader_drop_description_support(self.el)

    def require_multiseq_support(self):
        self.ssptab = True
        gtlib.gt_encseq_loader_require_multiseq_support(self.el)

    def drop_multiseq_support(self):
        self.ssptab = False
        gtlib.gt_encseq_loader_drop_multiseq_support(self.el)

    def require_des_tab(self):
        self.destab = True
        gtlib.gt_encseq_loader_require_des_tab(self.el)

    def do_not_require_des_tab(self):
        self.destab = False
        gtlib.gt_encseq_loader_do_not_require_des_tab(self.el)

    def require_ssp_tab(self):
        self.ssptab = True
        gtlib.gt_encseq_loader_require_ssp_tab(self.el)

    def do_not_require_ssp_tab(self):
        self.ssptab = False
        gtlib.gt_encseq_encoder_do_not_require_ssp_tab(self.el)

    def require_sds_tab(self):
        self.sdstab = True
        gtlib.gt_encseq_loader_require_sds_tab(self.el)

    def do_not_require_sds_tab(self):
        self.sdstab = False
        gtlib.gt_encseq_loader_do_not_require_sds_tab(self.el)

    def load(self, indexname):
        if not os.path.exists(indexname + ".esq"):
            raise IOError("file not found: %s" % indexname + ".esq")
        if self.destab and not os.path.exists(indexname + ".des"):
            raise IOError("file not found: %s" % indexname + ".des")
        # not required in every case (equallength seqs)
        # if self.ssptab and not os.path.exists(indexname+".ssp"):
        #    raise IOError, ("file not found: %s" % indexname+".ssp")
        if self.sdstab and not os.path.exists(indexname + ".sds"):
            raise IOError("file not found: %s" % indexname + ".sds")
        err = Error()
        esptr = gtlib.gt_encseq_loader_load(self.el,
                                            str(indexname).encode("UTF-8"),
                                            err.error)
        if not esptr:
            gterror(err)
        return Encseq(esptr, True)

    def register(cls, gtlib):
        gtlib.gt_encseq_loader_new.restype = c_void_p
        gtlib.gt_encseq_loader_load.argtypes = [c_void_p, c_char_p,  c_void_p]
        gtlib.gt_encseq_loader_load.restype = c_void_p
        gtlib.gt_encseq_loader_delete.argtypes = [c_void_p]

    register = classmethod(register)


class EncseqBuilder:

    def __init__(self, a):
        self.eb = gtlib.gt_encseq_builder_new(a.alpha)
        self._as_parameter_ = self.eb

    def __del__(self):
        try:
            gtlib.gt_encseq_builder_delete(self.eb)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, EncseqBuilder):
            raise TypeError("argument must be an EncseqBuilder object")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def enable_description_support(self):
        gtlib.gt_encseq_builder_enable_description_support(self.eb)

    def disable_description_support(self):
        gtlib.gt_encseq_builder_disable_description_support(self.eb)

    def enable_multiseq_support(self):
        gtlib.gt_encseq_builder_enable_multiseq_support(self.eb)

    def disable_multiseq_support(self):
        gtlib.gt_encseq_builder_disable_multiseq_support(self.eb)

    def create_des_tab(self):
        gtlib.gt_encseq_builder_create_des_tab(self.eb)

    def do_not_create_des_tab(self):
        gtlib.gt_encseq_builder_do_not_create_des_tab(self.eb)

    def create_ssp_tab(self):
        gtlib.gt_encseq_builder_create_ssp_tab(self.eb)

    def do_not_create_ssp_tab(self):
        gtlib.gt_encseq_builder_do_not_create_ssp_tab(self.eb)

    def create_sds_tab(self):
        gtlib.gt_encseq_builder_create_sds_tab(self.eb)

    def do_not_create_sds_tab(self):
        gtlib.gt_encseq_builder_do_not_create_sds_tab(self.eb)

    def add_string(self, string, desc=''):
        string = str(string)
        gtlib.gt_encseq_builder_add_cstr(self.eb, string.encode('UTF-8'),
                                         len(string), desc.encode('UTF-8'))

    def build(self):
        err = Error()
        esptr = gtlib.gt_encseq_builder_build(self.eb, err.error)
        if not esptr:
            gterror(err)
        return Encseq(esptr, True)

    def register(cls, gtlib):
        gtlib.gt_encseq_builder_new.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_new.restype = c_void_p
        gtlib.gt_encseq_builder_add_cstr.argtypes = [c_void_p, c_char_p,
                                                     c_ulong, c_char_p]
        gtlib.gt_encseq_builder_enable_description_support.argtypes = [
            c_void_p]
        gtlib.gt_encseq_builder_disable_description_support.argtypes = [
            c_void_p]
        gtlib.gt_encseq_builder_enable_multiseq_support.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_disable_multiseq_support.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_create_des_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_do_not_create_des_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_create_ssp_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_do_not_create_ssp_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_create_sds_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_do_not_create_sds_tab.argtypes = [c_void_p]
        gtlib.gt_encseq_builder_build.argtypes = [c_void_p, c_void_p]
        gtlib.gt_encseq_builder_build.restype = c_void_p
        gtlib.gt_encseq_builder_delete.argtypes = [c_void_p]

    register = classmethod(register)


class EncseqReader:

    def __init__(self, ptr, own=False):
        self.er = ptr
        self._as_parameter_ = self.er
        self.own = own

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_encseq_reader_delete(self.er)
            except AttributeError:
                pass

    def from_param(cls, obj):
        if not isinstance(obj, EncseqReader):
            raise TypeError("argument must be an EncseqReader")
        return obj._as_parameter_

    def next_encoded_char(self):
        return gtlib.gt_encseq_reader_next_encoded_char(self.er)

    def next_decoded_char(self):
        return gtlib.gt_encseq_reader_next_decoded_char(self.er)

    from_param = classmethod(from_param)

    def register(cls, gtlib):
        gtlib.gt_encseq_reader_next_encoded_char.restype = c_ubyte
        gtlib.gt_encseq_reader_next_encoded_char.argtypes = [c_void_p]
        gtlib.gt_encseq_reader_next_decoded_char.restype = c_char
        gtlib.gt_encseq_reader_next_decoded_char.argtypes = [c_void_p]
        gtlib.gt_encseq_reader_delete.argtypes = [c_void_p]

    register = classmethod(register)


class Encseq:

    def __init__(self, ptr, own=False):
        self.encseq = ptr
        self._as_parameter_ = self.encseq
        self.own = own

    def __del__(self):
        if self.own:
            try:
                gtlib.gt_encseq_delete(self.encseq)
            except AttributeError:
                pass

    def from_param(cls, obj):
        if not isinstance(obj, Encseq):
            raise TypeError("argument must be an Encseq")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def num_of_sequences(self):
        return gtlib.gt_encseq_num_of_sequences(self.encseq)

    def num_of_files(self):
        return gtlib.gt_encseq_num_of_files(self.encseq)

    def description(self, num):
        if not num < self.num_of_sequences():
            gterror("invalid sequence number %d" % num)
        desclen = c_ulong()
        str_addr = gtlib.gt_encseq_description(
            self.encseq, byref(desclen), num)
        return str(string_at(str_addr, desclen.value).decode('UTF-8'))

    def total_length(self):
        return gtlib.gt_encseq_total_length(self.encseq)

    def get_encoded_char(self, pos, readmode=0):
        if readmode < 0 or readmode > 3:
            gterror("invalid readmode!")
        return gtlib.gt_encseq_get_encoded_char(self.encseq, pos, readmode)

    def get_decoded_char(self, pos, readmode=0):
        if readmode < 0 or readmode > 3:
            gterror("invalid readmode!")
        return gtlib.gt_encseq_get_decoded_char(self.encseq, pos, readmode)

    def alphabet(self):
        a = gtlib.gt_encseq_alphabet(self.encseq)
        if not a:
            return None
        return Alphabet(a)

    def seqstartpos(self, num):
        if not num < self.num_of_sequences():
            gterror("invalid sequence number %d" % num)
        return gtlib.gt_encseq_seqstartpos(self.encseq, num)

    def filestartpos(self, num):
        if not num < self.num_of_files():
            gterror("invalid file number %d" % num)
        return gtlib.gt_encseq_filestartpos(self.encseq, num)

    def seqlength(self, num):
        if not num < self.num_of_sequences():
            gterror("invalid sequence number %d" % num)
        return gtlib.gt_encseq_seqlength(self.encseq, num)

    def effective_filelength(self, num):
        if not num < self.num_of_files():
            gterror("invalid file number %d" % num)
        return gtlib.gt_encseq_effective_filelength(self.encseq, num)

    def filenames(self, num):
        arr = StrArray(gtlib.gt_encseq_filenames(seld.encseq))
        return arr.to_list()

    def create_reader_with_readmode(self, readmode, startpos):
        if readmode < 0 or readmode > 3:
            gterror("invalid readmode!")
        if startpos < 0 or startpos >= self.total_length():
            gterror("invalid startposition: %d (allowed: %d-%d)" % (startpos,
                                                                    0,
                                                                    self.total_length() - 1))
        er = gtlib.gt_encseq_create_reader_with_readmode(self.encseq, readmode,
                                                         startpos)
        return EncseqReader(er, True)

    def extract_encoded(self, start, end):
        if start < 0 or end >= self.total_length():
            gterror("invalid coordinates: %d-%d (allowed: %d-%d)" % (start, end,
                                                                     0, self.total_length() - 1))
        buf = (c_ubyte * (end - start + 1))()
        gtlib.gt_encseq_extract_encoded(self.encseq, buf, start, end)
        return buf

    def extract_decoded(self, start, end):
        if start < 0 or end >= self.total_length():
            gterror("invalid coordinates: %d-%d (allowed: %d-%d)" % (start, end,
                                                                     0, self.total_length() - 1))
        buf = (c_char * (end - start + 1))()
        gtlib.gt_encseq_extract_decoded(self.encseq, buf, start, end)
        return string_at(buf, end - start + 1).decode('UTF-8')

    def seqnum(self, pos):
        if pos < 0 or pos >= self.total_length():
            gterror("invalid position: %d (allowed: %d-%d)" % (pos,
                                                               0, self.total_length() - 1))
        return gtlib.gt_encseq_seqnum(self.encseq, pos)

    def filenum(self, pos):
        if pos < 0 or pos >= self.total_length():
            gterror("invalid position: %d (allowed: %d-%d)" % (pos,
                                                               0, self.total_length() - 1))
        return gtlib.gt_encseq_filenum(self.encseq, pos)

    def mirror(self):
        if self.is_mirrored():
            gterror("encoded sequence is already mirrored")
        err = Error()
        ret = gtlib.gt_encseq_mirror(self.encseq, err.error)
        if ret != 0:
            gterror(err)

    def unmirror(self):
        if not self.is_mirrored():
            gterror("encoded sequence is not mirrored")
        err = Error()
        gtlib.gt_encseq_unmirror(self.encseq)

    def is_mirrored(self):
        return int2bool(gtlib.gt_encseq_is_mirrored(self.encseq))

    def register(cls, gtlib):
        gtlib.gt_encseq_create_reader_with_readmode.restype = c_void_p
        gtlib.gt_encseq_create_reader_with_readmode.argtypes = [
            c_void_p, c_int, c_ulong]
        gtlib.gt_encseq_num_of_sequences.restype = c_ulong
        gtlib.gt_encseq_num_of_sequences.argtypes = [c_void_p]
        gtlib.gt_encseq_num_of_files.restype = c_ulong
        gtlib.gt_encseq_num_of_files.argtypes = [c_void_p]
        gtlib.gt_encseq_total_length.restype = c_ulong
        gtlib.gt_encseq_total_length.argtypes = [c_void_p]
        gtlib.gt_encseq_description.restype = c_char_p
        gtlib.gt_encseq_description.argtypes = [c_void_p, POINTER(c_ulong),
                                                c_ulong]
        gtlib.gt_encseq_get_encoded_char.restype = c_ubyte
        gtlib.gt_encseq_get_encoded_char.argtypes = [c_void_p, c_ulong, c_int]
        gtlib.gt_encseq_get_decoded_char.restype = c_char
        gtlib.gt_encseq_get_decoded_char.argtypes = [c_void_p, c_ulong, c_int]
        gtlib.gt_encseq_alphabet.restype = c_void_p
        gtlib.gt_encseq_alphabet.argtypes = [c_void_p]
        gtlib.gt_encseq_seqstartpos.restype = c_ulong
        gtlib.gt_encseq_seqstartpos.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_seqlength.restype = c_ulong
        gtlib.gt_encseq_seqlength.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_effective_filelength.restype = c_uint64
        gtlib.gt_encseq_effective_filelength.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_extract_encoded.argtypes = [c_void_p,
                                                    POINTER(c_ubyte),
                                                    c_ulong, c_ulong]
        gtlib.gt_encseq_extract_decoded.argtypes = [c_void_p, c_char_p,
                                                    c_ulong, c_ulong]
        gtlib.gt_encseq_mirror.restype = c_int
        gtlib.gt_encseq_mirror.argtypes = [c_void_p, c_void_p]
        gtlib.gt_encseq_unmirror.argtypes = [c_void_p]
        gtlib.gt_encseq_is_mirrored.restype = c_int
        gtlib.gt_encseq_is_mirrored.argtypes = [c_void_p]
        gtlib.gt_encseq_filestartpos.restype = c_ulong
        gtlib.gt_encseq_filestartpos.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_seqnum.restype = c_ulong
        gtlib.gt_encseq_seqnum.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_filenum.restype = c_ulong
        gtlib.gt_encseq_filenum.argtypes = [c_void_p, c_ulong]
        gtlib.gt_encseq_filenames.restype = c_void_p
        gtlib.gt_encseq_filenames.argtypes = [c_void_p]
        gtlib.gt_encseq_delete.argtypes = [c_void_p]

    register = classmethod(register)
