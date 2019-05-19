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

require 'dl/import'
require 'dl/struct'
require 'gthelper'
require 'extended/feature_node'
require 'annotationsketch/rec_map'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtImageInfo* gt_image_info_new()"
  extern "unsigned int gt_image_info_get_height(GtImageInfo*)"
  extern "unsigned long gt_image_info_num_of_rec_maps(GtImageInfo*)"
  extern "const GtRecMap* gt_image_info_get_rec_map(GtImageInfo*, " +
                                                   "unsigned long)"
  extern "void gt_image_info_delete(GtImageInfo*)"

  class ImageInfo
    attr_reader :image_info
    def initialize()
      @image_info = GT.gt_image_info_new()
      @image_info.free = GT::symbol("gt_image_info_delete", "0P")
    end

    def get_height()
      GT.gt_image_info_get_height(@image_info)
    end

    def num_of_rec_maps()
      GT.gt_image_info_num_of_rec_maps(@image_info)
    end

    def each_hotspot()
      if @hotspots.nil? then
        @hotspots = []
        0.upto(self.num_of_rec_maps()-1) do |i|
          rm = GT::RecMap.new(GT.gt_image_info_get_rec_map(@image_info, i))
          @hotspots.push([rm.get_northwest_x.to_i, \
                          rm.get_northwest_y.to_i, \
                          rm.get_southeast_x.to_i, \
                          rm.get_southeast_y.to_i, \
                          rm.get_genome_feature])
        end
        @hotspots.sort!{|hs1,hs2|
          if hs1[2]-hs1[0]+1 == hs2[2]-hs2[0]+1 then
            hs1[3] <=> hs2[3]
          else
            hs1[2]-hs1[0]+1 <=> hs2[2]-hs2[0]+1
          end
        }
      end
      @hotspots.each do |hs|
        yield hs[0],hs[1],hs[2],hs[3],hs[4]
      end
    end

    def to_ptr
      @image_info
    end
  end
end
