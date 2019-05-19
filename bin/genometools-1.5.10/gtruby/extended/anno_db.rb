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
require 'extended/feature_index'

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  extern "GtAnnoDBSchema* gt_anno_db_gfflike_new()"
  extern "GtFeatureIndex* gt_anno_db_schema_get_feature_index(GtAnnoDBSchema*,
                                                              GtRDB*,
                                                              GtError*)"
  extern "void            gt_anno_db_schema_delete(GtAnnoDBSchema*) "

  class AnnoDBSchema
    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end

    def to_ptr
      @adb
    end

    def get_feature_index(rdb)
      err = GT::Error.new()
      fi = GT.gt_anno_db_schema_get_feature_index(@adb, rdb.to_ptr, err)
      if fi.nil? then
        GT::gterror(err)
      end
      GT::FeatureIndexFromPtr.new(fi)
    end
  end

  class AnnoDBGFFlike < AnnoDBSchema
    def initialize()
      err = GT::Error.new()
      @adb = GT.gt_anno_db_gfflike_new()
      if @adb.nil? then
        GT::gterror(err)
      end
      @adb.free = GT::symbol("gt_anno_db_schema_delete", "0P")
    end
  end
end
