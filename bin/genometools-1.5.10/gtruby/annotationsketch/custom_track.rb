#
# Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg
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
  extern "GtCustomTrack* gt_custom_track_script_wrapper_new(void*, void*, void*,
                                                            void*)"

  class CustomTrack
    def initialize()
      @get_height = DL.callback("L") do
                      ret = self.get_height()
                      if ret.nil? or !ret.is_a?(Numeric) then
                        GT::gterror("custom track get_height() method must " + \
                                    "return a numeric value!")
                      end
                      ret
                    end
      @get_title = DL.callback("0PP") do |ptr, str_ptr|
                      ret = self.get_title()
                      s = GT::Str.new(str_ptr)
                      if ret.nil? or !ret.to_s.is_a?(String) then
                        GT::gterror("custom track get_title() method must " +  \
                                    "return a string!")
                      end
                      s.append_cstr(ret)
                    end
      @render    = DL.callback("IPLPPP") do |g, ypos, rng, sty, err|
                      rng.struct!("LL", :start, :end)
                      range = (rng[:start]..rng[:end])
                      ret = self.render(GT::Graphics.new(g), ypos, range,        \
                                        GT::Style.new(sty), GT::Error.new(err))
                      if ret.nil? or !ret.is_a?(Numeric) then
                        GT::gterror("custom track render() method must " +     \
                                    "return a numeric value!")
                      end
                      ret
                    end
      @free    = DL.callback("0") do
                      self.free()
                    end
      @ct = GT.gt_custom_track_script_wrapper_new(@render,                     \
                                                  @get_height,                 \
                                                  @get_title,                  \
                                                  @free)
      @ct.free = GT::symbol("gt_custom_track_delete", "0P")
    end

    def to_ptr
      @ct
    end
  end
end
