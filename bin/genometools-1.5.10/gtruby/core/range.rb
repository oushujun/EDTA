#
# Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
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

module GT
  class Range
    def initialize(start, stop)
      set(start, stop)
    end

    def set(start, stop)
      if start.nil? or stop.nil? or start > stop or start < 0 or stop < 0 then
        GT.gterror("range error: start > end!")
      end
      @start = start
      # 'end' is a reserved token in Ruby and cannot be used as an identifier
      @stop = stop
    end

    def start
      @start
    end

    def start=(val)
      if val > @stop or not val >= 0 then
        GT.gterror("Invalid range start component: %d" % val)
      end
      @start = val
    end

    def Range.malloc
      Range.new()
    end

    def end
      @stop
    end
    
    def end=(val)
      if val < @start or not val >= 0 then
        GT.gterror("Invalid range end component: %d" % val)
      end
      @stop = val
    end

    def to_ptr
      (@mem = [@start, @stop].pack("L_L_")).to_ptr
    end
  end
end
