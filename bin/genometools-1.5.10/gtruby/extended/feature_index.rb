#
# Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg
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
require 'core/array'
require 'core/range'
require 'core/str_array'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  extern "GtFeatureIndex* gt_feature_index_memory_new()"
  extern "void gt_feature_index_delete(GtFeatureIndex*)"
  extern "int gt_feature_index_add_gff3file(GtFeatureIndex*, " +
                                           "const char*, " +
                                           "GtError*)"
  extern "GtArray* gt_feature_index_get_features_for_seqid(Gt_FeatureIndex*, " +
                                                          "const char*, " +
                                                          "GtError*)"
  extern "char* gt_feature_index_get_first_seqid(const GtFeatureIndex*, " +
                                                "GtError*)"
  extern "GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex*, " +
                                                 "GtError*)"
  extern "int  gt_feature_index_get_range_for_seqid(GtFeatureIndex*, " +
                                                   "GtRange*, const char*, " +
                                                   "GtError*)"
  extern "int  gt_feature_index_has_seqid(const GtFeatureIndex*, int*, " +
                                        "const char*, GtError*)"
  extern "int  gt_feature_index_get_features_for_range(GtFeatureIndex*,
                                                       GtArray*, const char*,
                                                       const GtRange*,
                                                       GtError*)"
  class FeatureIndex
    attr_reader :feature_index

    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end

    def get_features_for_seqid(seqid)
      err = GT::Error.new()
      rval = GT.gt_feature_index_get_features_for_seqid(@feature_index, seqid, \
                                                        err)
      if rval then
        a = GT::Array.new(rval)
        result = []
        1.upto(a.size) do |i|
          fn = GT::FeatureNode.new(GT.gt_genome_node_ref(a.get(i-1)))
          result.push(fn)
        end
        result
      else
        GT.gterror(err)
      end
    end

    def get_features_for_range(start, stop, seqid)
      a = Array.create()
      err = Error.new()
      rng = Range.new(start, stop)
      rval = GT::gt_feature_index_get_features_for_range(@feature_index, a, \
                                                         seqid, rng.to_ptr, err)
      if rval != 0 then
          GT.gterror(err)
      end
      result = []
      1.upto(a.size) do |i|
        fn = GT::FeatureNode.new(GT.gt_genome_node_ref(a.get(i-1)))
        result.push(fn)
      end
      result
    end

    def add_gff3file(filename)
      err = GT::Error.new()
      rval = GT.gt_feature_index_add_gff3file(@feature_index, filename, \
                                              err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def get_first_seqid
      err = Error.new()
      val = GT.gt_feature_index_get_first_seqid(@feature_index, err.to_ptr)
      if val.nil? then
        GT.gterror(err) if err.is_set?
      end
      val
    end

    def get_seqids
      err = Error.new()
      sap = GT.gt_feature_index_get_seqids(@feature_index, err.to_ptr)
      if sap.nil? then GT.gterror(err) end
      GT::StrArray.new(sap).to_a
    end

    def has_seqid?(seqid)
      err = Error.new()
      val = (GT::struct ["bool val"]).malloc
      rval = GT.gt_feature_index_has_seqid(@feature_index, val, seqid, err)
      if rval != 0 then
        GT.gterror(err)
      else
        return (val.val == true)
      end
    end

    def get_range_for_seqid(seqid)
      err = Error.new()
      if not self.has_seqid?(seqid) then
        GT.gterror("feature_index does not contain seqid")
      end
      range = DL::malloc(DL::sizeof("LL"))
      rval = GT.gt_feature_index_get_range_for_seqid(@feature_index, range, \
                                                     seqid, err)
      if rval != 0 then
        GT.gterror(err)
      else
        range.struct!("LL", :start, :stop)
        Range.new(range[:start],range[:stop])
      end
    end

    def to_ptr
      @feature_index
    end
  end

  class FeatureIndexMemory < FeatureIndex
    def initialize
      @feature_index = GT.gt_feature_index_memory_new()
      @feature_index.free = GT::symbol("gt_feature_index_delete", "0P")
    end
  end

  class FeatureIndexFromPtr < FeatureIndex
    def initialize(ptr, own = false)
      @feature_index = ptr
      if own then
        @feature_index.free = GT::symbol("gt_feature_index_delete", "0P")
      end
    end
  end
end
