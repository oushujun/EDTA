#
# Copyright (c) 2009-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009-2012 Center for Bioinformatics, University of Hamburg
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
require 'extended/comment_node'
require 'extended/eof_node'
require 'extended/feature_node'
require 'extended/meta_node'
require 'extended/sequence_node'
require 'extended/region_node'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtNodeVisitor* gt_script_wrapper_visitor_new(void*, void*, void*,
                                                       void*, void*, void*,
                                                       void*)"

  class CustomVisitor
    def visitor_func_generic(node_class, node_ptr, err_ptr, method_name)
      node = node_class.new(node_ptr, true)
      err = GT::Error.new(err_ptr)
      begin
        begin
          ret = self.send("visit_#{method_name}_node", node)
        rescue NoMethodError
          ret = 0
        end
      rescue Error => msg
        err.set(msg)
        ret = -1
      end
      if ret.nil?
        0
      else
        begin
          Integer(ret)
        rescue
          err.set("Function 'visit_#{method_name}_node' of class " + \
                  "'#{self.class}' must return an integer or nil!")
          1
        end
      end
    end

    def initialize()
      # try to save callbacks by registering them only for implemented handlers
      if self.respond_to?("visit_feature_node") then
        @feature_node_cb = DL.callback("IPP") do |fn_ptr, err_ptr|
          self.visitor_func_generic(GT::FeatureNode, fn_ptr, err_ptr, "feature")
        end
      end
      if self.respond_to?("visit_comment_node") then
        @comment_node_cb = DL.callback("IPP") do |cn_ptr, err_ptr|
          self.visitor_func_generic(GT::CommentNode, cn_ptr, err_ptr, "comment")
        end
      end
      if self.respond_to?("visit_region_node") then
        @region_node_cb = DL.callback("IPP") do |rn_ptr, err_ptr|
          self.visitor_func_generic(GT::RegionNode, rn_ptr, err_ptr, "region")
        end
      end
      if self.respond_to?("visit_sequence_node") then
        @sequence_node_cb = DL.callback("IPP") do |sn_ptr, err_ptr|
          self.visitor_func_generic(GT::SequenceNode, sn_ptr, err_ptr, "sequence")
        end
      end
      if self.respond_to?("visit_meta_node") then
        @meta_node_cb = DL.callback("IPP") do |mn_ptr, err_ptr|
          self.visitor_func_generic(GT::MetaNode, mn_ptr, err_ptr, "meta")
        end
      end
      if self.respond_to?("visit_eof_node") then
        @eof_node_cb = DL.callback("IPP") do |en_ptr, err_ptr|
          self.visitor_func_generic(GT::EOFNode, en_ptr, err_ptr, "eof")
        end
      end
      @genome_visitor = GT.gt_script_wrapper_visitor_new(@comment_node_cb,
                                                         @feature_node_cb,
                                                         @region_node_cb,
                                                         @sequence_node_cb,
                                                         @meta_node_cb,
                                                         @eof_node_cb,
                                                         nil)
      @genome_visitor.free = GT::symbol("gt_node_visitor_delete", "0P")
    end

    def to_ptr
      @genome_visitor
    end
  end
end
