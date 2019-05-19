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

module GT
  extend DL::Importable
  gtdlload "libgenometools"

  extern "GtRDB* gt_rdb_sqlite_new(const char*, GtError*)"
  begin
    extern "GtRDB* gt_rdb_mysql_new(const char*, unsigned int,
                                    const char*, const char*,
                                    const char*, GtError*)"
  rescue
  end

  class RDB
    def initialize(*)
      raise(NotImplementedError, "Please call the constructor of a " +
                                 "#{self.class} implementation.")
    end

    def to_ptr
      @rdb
    end
  end

  class RDBSqlite < RDB
    def initialize(dbpath)
      err = GT::Error.new()
      @rdb = GT.gt_rdb_sqlite_new(dbpath, err)
      if @rdb.nil? then
        GT::gterror(err)
      end
      @rdb.free = GT::symbol("gt_rdb_delete", "0P")
    end
  end

  class RDBMySQL < RDB
    def initialize(server, port, db, user, pass)
      err = GT::Error.new()
      @rdb = GT.gt_rdb_mysql_new(server, port, db, user, pass, err)
      if @rdb.nil? then
        GT::gterror(err)
      end
      @rdb.free = GT::symbol("gt_rdb_delete", "0P")
    end
  end
end
