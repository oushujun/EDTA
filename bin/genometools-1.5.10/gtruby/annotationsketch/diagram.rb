#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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
require 'annotationsketch/block'
require 'core/array'
require 'core/range'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtDiagram* gt_diagram_new(GtFeatureIndex*, const char*, " + \
                                    "const GtRange*, GtStyle*, GtError*)"
  extern "GtDiagram* gt_diagram_new_from_array(GtArray*, " + \
                                               "const GtRange*, GtStyle*)"
  extern "void gt_diagram_set_track_selector_func(GtDiagram*, void*, void*)"
  extern "void gt_diagram_add_custom_track(GtDiagram*, GtCustomTrack*)"
  extern "void gt_diagram_delete(GtDiagram*)"

  class Diagram
    attr_reader :diagram

    def self.from_index(feature_index, seqid, range, style)
      err = GT::Error.new()
      if range.start > range.end
        GT.gterror("range.start > range.end")
      end
      if feature_index.nil? then
        GT.gterror("feature index must not be nil!")
      end
      if seqid.nil? then
        GT.gterror("seqid must not be nil!")
      end
      if !style.is_a?(GT::Style) then
        GT.gterror("'style' parameter must be a Style object!")
      end
      diagram = GT.gt_diagram_new(feature_index, seqid, range.to_ptr,
                                  style, err)
      if diagram.nil? then
        GT::gterror(err)
      end
      return GT::Diagram.new(diagram)
    end

    def self.from_array(array, range, style)
      if range.start > range.end
        GT.gterror("range.start > range.end")
      end
      gtarr = GT::Array.create(DL.sizeof("P"), false)
      array.each do |i|
        if !i.is_a?(GT::FeatureNode) then
          GT::gterror("Diagram array must only contain FeatureNodes!")
        end
        gtarr.add(i)
      end
      diagram = GT.gt_diagram_new_from_array(gtarr, range.to_ptr, style)
      return GT::Diagram.new(diagram)
    end

    def initialize(ptr)
      @diagram = ptr
      @diagram.free = GT::symbol("gt_diagram_delete", "0P")
#      callback_releaser = Proc.new do
#        DL.remove_callback(@tsf) unless @tsf.nil
#      end
#      ObjectSpace.define_finalizer(self, callback_releaser)
    end

    def set_track_selector_func(proc)
      @tsf = DL.callback('0PPP') do |b_ptr, string_ptr, data_ptr|
               b = GT::Block.new(b_ptr)
               string = GT::Str.new(string_ptr)
               s = proc.call(b, data_ptr)
               if (s.nil?) then
                 GT::gterror("Track selector callback must return a string!")
               else
                 string.append_cstr(s)
               end
             end
      GT.gt_diagram_set_track_selector_func(@diagram, @tsf, GT::NULL)
    end

    def release_track_selector_func
      DL.remove_callback(@tsf)
    end

    def add_custom_track(ct)
      GT.gt_diagram_add_custom_track(@diagram, ct)
    end

    def to_ptr
      @diagram
    end
  end
end
