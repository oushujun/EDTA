#
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
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

require 'gtdlload'
require 'core/str'
require 'core/str_array'

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  typealias "bool", "ibool"
  typealias "GtUchar", "unsigned char"
  extern "GtAlphabet* gt_alphabet_new_from_sequence(const GtStrArray*,
                                                    GtError*)"
  extern "GtAlphabet* gt_alphabet_new_dna()"
  extern "GtAlphabet* gt_alphabet_new_protein()"
  extern "GtAlphabet* gt_alphabet_new_empty()"
  extern "GtAlphabet* gt_alphabet_new_from_file(const char*, GtError*)"
  extern "GtAlphabet* gt_alphabet_guess(const char*, unsigned long)"
  extern "GtAlphabet* gt_alphabet_ref(GtAlphabet*)"
  extern "void gt_alphabet_delete(GtAlphabet*)"
  extern "unsigned int gt_alphabet_num_of_chars(const GtAlphabet*)"
  extern "unsigned int gt_alphabet_size(const GtAlphabet*)"
  extern "const void* gt_alphabet_characters(const GtAlphabet*)"
  extern "int gt_alphabet_is_protein(const GtAlphabet*)"
  extern "int gt_alphabet_is_dna(const GtAlphabet*)"
  extern "int gt_alphabet_valid_input(const GtAlphabet*, char)"
  extern "GtUchar gt_alphabet_encode(const GtAlphabet*, char)"
  extern "char gt_alphabet_decode(const GtAlphabet*, GtUchar)"
  extern "GtStr* gt_alphabet_decode_seq_to_str(const GtAlphabet*,
                                               const GtUchar*,
                                               unsigned long)"

  class Alphabet

    def self.create_dna
      return Alphabet.new(GT.gt_alphabet_new_dna())
    end

    def self.create_protein
      return Alphabet.new(GT.gt_alphabet_new_protein())
    end

    def self.create_from_file(filename)
      err = Error.new
      aptr = GT.gt_alphabet_new_from_file(filename.to_str, err.to_ptr)
      if aptr == GT::NULL
        GT.gterror(err)
      else
        return Alphabet.new(aptr)
      end
    end

    def self.create_from_sequence(files)
      files = files.to_a
      if files.length == 0 then
        GT.gterror("list of input sequence files must be non-empty")
      end
      sa = StrArray.new
      files.each do |f|
        fn = f.to_s
        if (!File.exists?(fn))
          GT.gterror("file not found: #{fn}")
        end
        sa.add(fn)
      end
      err = Error.new
      aptr = GT.gt_alphabet_new_from_sequence(sa, err)
      if aptr == GT::NULL
        GT.gterror(err)
      else
        return Alphabet.new(aptr)
      end
    end

    def initialize(aptr, own = true)
      if own then
        @alpha = aptr
      else
        @alpha = GT.gt_alphabet_ref(aptr)
      end
      @alpha.free = GT::symbol("gt_alphabet_delete", "0P")
    end

    def size
      GT.gt_alphabet_size(@alpha)
    end

    def num_of_chars
      GT.gt_alphabet_num_of_chars(@alpha)
    end

    def valid_input(char)
      GT.gt_alphabet_valid_input(@alpha, char) != 0
    end

    def is_dna?
      GT.gt_alphabet_is_dna(@alpha) != 0
    end

    def is_protein?
      GT.gt_alphabet_is_protein(@alpha) != 0
    end

    def characters
      ptr = GT.gt_alphabet_characters(@alpha)
      ptr.to_a('C', self.num_of_chars).collect{|c| c.chr}
    end

    def encode(char)
      if !self.valid_input(char) then
        GT.gterror("#{char} is not a valid input character!")
      end
      return GT.gt_alphabet_encode(@alpha, char)
    end

    def decode(char)
      return GT.gt_alphabet_decode(@alpha, char).chr
    end

    def decode_seq(arr)
      strptr = GT.gt_alphabet_decode_seq_to_str(@alpha, arr.pack('C*'), \
                                                arr.length)
      Str.new(strptr).to_s
    end

    def to_ptr
      @alpha
    end
  end
end
