#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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

require 'dl/import'
require 'gthelper'
require "core/str"

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtGenomeNode* gt_sequence_node_new(const char*, GtStr*)"
  extern "const char*   gt_sequence_node_get_description(const GtSequenceNode*)"
  extern "const char*   gt_sequence_node_get_sequence(const GtSequenceNode*)"
  extern "unsigned long gt_sequence_node_get_sequence_length(const
                                                             GtSequenceNode*)"

  class SequenceNode < GenomeNode
    def self.create(description, sequence)
      newsn = GT.gt_sequence_node_new(description.to_s, \
                                      GT::Str.new(sequence.to_s))
      return GT::SequenceNode.new(newsn, true)
    end

    def initialize(gn, newref=false)
      super(gn, newref)
    end

    def get_description
      return GT.gt_sequence_node_get_description(@genome_node)
    end

    def get_sequence
      (seq = GT.gt_sequence_node_get_sequence(@genome_node)) ? seq : ""
    end

    def get_sequence_length
      return GT.gt_sequence_node_get_sequence_length(@genome_node)
    end
  end
end
