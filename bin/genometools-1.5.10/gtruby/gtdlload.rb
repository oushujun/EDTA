#
# Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

module GT
  extend DL::Importable

  # XXX: Does anybody have a better idea for this mess?
  def GT.gtdlload(basename)
    if not $GT_SYSTEM then
      $GT_SYSTEM=`uname -s`.chomp
    end
    if $GT_SYSTEM == "Darwin" then
      dlload basename + ".dylib"
    else
      dlload basename + ".so"
    end
  end

  # a NULL pointer
  NULL = DL::PtrData.new(0)
end

if __FILE__ == $0
  GT.gtdlload "libgenometools"
end
