#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg
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
require 'gthelper'
require 'core/error'
require 'extended/genome_node'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "int gt_node_stream_next(GtNodeStream*, GtGenomeNode**, GtError*)"

  class GenomeStream

    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end

    def next_tree
      err = GT::Error.new()
      genome_node = DL::PtrData.new(0)
      genome_node.free = DL::FREE
      rval = GT.gt_node_stream_next(@genome_stream, genome_node.ref,
                                    err.to_ptr)
      if rval != 0 then GT.gterror(err) end
      if genome_node.null? then return nil end
      GT::GenomeNode.new(genome_node)
    end

    def to_ptr
      @genome_stream
    end
  end
end
