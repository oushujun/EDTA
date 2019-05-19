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

require 'gtdlload'
require 'gthelper'
require 'extended/genome_stream'
require 'dl/struct'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtNodeStream* gt_script_wrapper_stream_new(void*, void*)"

  class CustomStream < GenomeStream
    def initialize()
      if !public_methods.include?("next") then
        GT::gterror "#{self.class} must implement 'next' method"
      end

      @next_cb = DL.callback("IPP") do |gn_ptr, err_ptr|
        rval = 0
        begin
          node = self.next
          gn_ptr.struct!("P", :val)
          if !node.nil? then
            GT.gt_genome_node_ref(node)
            gn_ptr[:val] = node.to_ptr
          else
            gn_ptr[:val] = nil
          end
        rescue => errmsg
          Error.new(err_ptr).set(errmsg)
          rval = -1
        end
        rval
      end

      @genome_stream = GT.gt_script_wrapper_stream_new(@next_cb, nil)
      @genome_stream.free = GT::symbol("gt_node_stream_delete", "0P")
    end

    def to_ptr
      @genome_stream
    end
  end
end
