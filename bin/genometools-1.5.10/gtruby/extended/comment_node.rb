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

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtCommentNode* gt_comment_node_new(const char*)"
  extern "const char*    gt_comment_node_get_comment(const GtCommentNode*)"

  class CommentNode < GenomeNode
    def self.create(comment)
      newfn = GT.gt_comment_node_new(comment.to_s)
      return GT::CommentNode.new(newfn, true)
    end

    def initialize(gn, newref=false)
      super(gn, newref)
    end

    def get_comment
      return GT.gt_comment_node_get_comment(@genome_node)
    end
  end
end
