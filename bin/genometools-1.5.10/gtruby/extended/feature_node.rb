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
require 'gthelper'
require 'core/str'
require 'extended/genome_node'
require 'extended/strand'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  typealias "GtStrand", "int"
  typealias "GtPhase", "int"
  extern "GtGenomeNode* gt_feature_node_new(GtStr*, const char*,
                                            unsigned long, unsigned long,
                                            GtStrand)"
  extern "void          gt_feature_node_add_child(GtFeatureNode*,
                                                  GtFeatureNode*)"
  extern "const char*   gt_feature_node_get_source(GtFeatureNode*)"
  extern "void          gt_feature_node_set_source(GtFeatureNode*, GtStr*)"
  extern "const char*   gt_feature_node_get_type(GtFeatureNode*)"
  extern "bool          gt_feature_node_has_type(GtFeatureNode*, const char*)"
  extern "bool          gt_feature_node_score_is_defined(const GtFeatureNode*)"
  extern "void          gt_feature_node_get_score_p(const GtFeatureNode*,
                                                    void*)"
  extern "void          gt_feature_node_set_score_p(GtFeatureNode*, float*)"
  extern "void          gt_feature_node_unset_score(GtFeatureNode*)"
  extern "GtStrand      gt_feature_node_get_strand(GtFeatureNode*)"
  extern "void          gt_feature_node_set_strand(GtFeatureNode*, GtStrand)"
  extern "GtPhase       gt_feature_node_get_phase(GtFeatureNode*)"
  extern "void          gt_feature_node_set_phase(GtFeatureNode*, GtPhase)"
  extern "const char*   gt_feature_node_get_attribute(GtFeatureNode*,
                                                      const char*)"
  extern "void          gt_feature_node_add_attribute(GtFeatureNode*,
                                                      const char*, const char*)"
  extern "void          gt_feature_node_set_attribute(GtFeatureNode*,
                                                      const char*, const char*)"
  extern "void          gt_feature_node_remove_attribute(GtFeatureNode*,
                                                         const char*)"
  extern "void          gt_feature_node_foreach_attribute(GtFeatureNode*,
                                                          void*, void*)"

  #callback to populate attribute list
  def collect_attrib(tag, val, data)
    GT.gt_str_array_add_cstr(data, tag)
    GT.gt_str_array_add_cstr(data, val)
  end
  COLLECTFUNC = callback "void collect_attrib(const char*, const char*, void*)"

  class FeatureNode < GenomeNode
    def self.create(seqid, type, start, stop, strand)
      s = GT::Str.new(seqid)
      if !GT::STRANDCHARS.include?(strand)
        GT::gterror("Invalid strand: '#{strand}'")
      end
      newfn = GT.gt_feature_node_new(s, type, start, stop, \
                                     GT::STRANDCHARS.index(strand))
      return GT::FeatureNode.new(newfn, true)
    end

    def initialize(gn, newref=false)
      super(gn, newref)
    end

    def add_child(node)
      if self.get_seqid != node.get_seqid
        GT::gterror("nodes must have identical sequence regions! " + \
                    "(was: '#{self.get_seqid}' vs. '#{node.get_seqid}'")
      end
      GT.gt_feature_node_add_child(@genome_node, node)
    end

    def get_source
      GT.gt_feature_node_get_source(@genome_node)
    end

    def set_source(source)
      s = GT::Str.new(source.to_s)
      GT.gt_feature_node_set_source(@genome_node, s)
    end

    def update_attribs
      attribs = GT::StrArray.new
      GT.gt_feature_node_foreach_attribute(@genome_node, COLLECTFUNC, attribs)
      attr_a = attribs.to_a
      @attribs = {}
      while not attr_a.empty? do
        @attribs[attr_a.shift] = attr_a.shift
      end
    end

    def get_type
      GT.gt_feature_node_get_type(@genome_node)
    end

    def has_type?(type)
      GT.gt_feature_node_has_type(@genome_node, type.to_s)
    end

    def score_is_defined?
      GT.gt_feature_node_score_is_defined(@genome_node)
    end

    def get_score
      if GT.gt_feature_node_score_is_defined(@genome_node) then
        n = DL::malloc(GT::NATIVEFLTSIZE)
        GT.gt_feature_node_get_score_p(@genome_node, n)
        return n[0, n.size].unpack("F")[0]
      else
        nil
      end
    end

    def set_score(score)
      GT.gt_feature_node_set_score_p(@genome_node, [score.to_f].pack("F").to_ptr)
    end

    def unset_score
      GT.gt_feature_node_unset_score(@genome_node)
    end

    def get_strand
      GT::STRANDCHARS[GT.gt_feature_node_get_strand(@genome_node)]
    end

    def set_strand
      if !GT::STRANDCHARS.include?(strand)
        GT::gterror("Invalid strand: '#{strand}'")
      end
      GT.gt_feature_node_set_strand(GT::STRANDCHARS.index(strand))
    end

    def get_phase
      GT.gt_feature_node_get_phase(@genome_node)
    end

    def set_phase(phase)
      if phase.to_i > 3 then
        GT::gterror("Invalid phase: '#{phase}'")
      end
      GT.gt_feature_node_set_phase(@genome_node, phase.to_i)
    end

    def get_attribute(attrib)
      if @attribs.nil?
        self.update_attribs
      end
      @attribs[attrib]
    end

    def add_attribute(tag, val)
      if tag.to_s == "" or val.to_s == "" then
        GT::gterror("Attribute keys or values must not be empty!")
      end
      GT.gt_feature_node_add_attribute(@genome_node, tag.to_s, val.to_s)
      self.update_attribs
    end

    def set_attribute(tag, val)
      if tag.to_s == "" or val.to_s == "" then
        GT::gterror("Attribute keys or values must not be empty!")
      end
      GT.gt_feature_node_set_attribute(@genome_node, tag.to_s, val.to_s)
      self.update_attribs
    end

    def remove_attribute(tag)
      if tag.to_s == "" then
        GT::gterror("Attribute key must not be empty!")
      end
      GT.gt_feature_node_remove_attribute(@genome_node, tag.to_s)
      self.update_attribs
    end

    def each_attribute
      if @attribs.nil?
        self.update_attribs
      end
      @attribs.each_pair do |tag, val|
        yield tag, val
      end
    end

    def traverse(it)
      tfn = it.next
      while !tfn.nil? do
        yield tfn
        tfn = it.next
      end
    end

    def traverse_dfs
      self.traverse(GT::FeatureNodeIteratorDepthFirst.new(self)) do |n|
        yield n
      end
    end

    def traverse_direct
      self.traverse(GT::FeatureNodeIteratorDirect.new(self)) do |n|
        yield n
      end
    end

    def to_ptr
      @genome_node
    end
  end

  extern "GtFeatureNode* gt_feature_node_iterator_next(GtFeatureNodeIterator*)"
  extern "void           gt_feature_node_iterator_delete(GtFeatureNodeIterator*)"
  extern "GtFeatureNodeIterator* gt_feature_node_iterator_new(const
                                                              GtFeatureNode*)"
  extern "GtFeatureNodeIterator* gt_feature_node_iterator_new_direct(const
                                                              GtFeatureNode*)"
  class FeatureNodeIterator
    def next
      ret = GT.gt_feature_node_iterator_next(@i)
      if !ret.nil? then
        GT::FeatureNode.new(ret, true)
      else
        nil
      end
    end
  end

  class FeatureNodeIteratorDepthFirst < FeatureNodeIterator
    def initialize(node)
      @i = GT.gt_feature_node_iterator_new(node)
      @i.free = GT::symbol("gt_feature_node_iterator_delete", "0P")
    end
  end

  class FeatureNodeIteratorDirect < FeatureNodeIterator
    def initialize(node)
      @i = GT.gt_feature_node_iterator_new_direct(node)
      @i.free = GT::symbol("gt_feature_node_iterator_delete", "0P")
    end
  end
end
