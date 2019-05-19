#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
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

module GT
  extend DL::Importable
  typealias "bool", "ibool"

  extern "double gt_rec_map_get_northwest_x(const GtRecMap*)"
  extern "double gt_rec_map_get_northwest_y(const GtRecMap*)"
  extern "double gt_rec_map_get_southeast_x(const GtRecMap*)"
  extern "double gt_rec_map_get_southeast_y(const GtRecMap*)"
  extern "const GtFeatureNode* gt_rec_map_get_genome_feature(const GtRecMap*)"
  extern "bool gt_rec_map_has_omitted_children(const GtRecMap*)"

  class RecMap
    def initialize(rm)
      @rm = rm
    end

    def get_northwest_x
      GT::gt_rec_map_get_northwest_x(@rm)
    end

    def get_northwest_y
      GT::gt_rec_map_get_northwest_y(@rm)
    end

    def get_southeast_x
      GT::gt_rec_map_get_southeast_x(@rm)
    end

    def get_southeast_y
      GT::gt_rec_map_get_southeast_y(@rm)
    end

    def get_genome_feature
      #refcount only this FeatureNode!
      GT::FeatureNode.new(GT::gt_rec_map_get_genome_feature(@rm), true)
    end

    def has_omitted_children
      GT::gt_rec_map_has_omitted_children(@gn)
    end
  end
end
