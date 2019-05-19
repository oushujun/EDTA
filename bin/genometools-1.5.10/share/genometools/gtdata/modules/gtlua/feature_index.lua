--[[
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

module(..., package.seeall)

require "gtlua.genome_features"

-- XXX: remove if statement if libannotationsketch is always compiled in
if GenomeTools_feature_index then
  -- Computes the coverage for the sequence ID <seqid>. The optional <maxdist>
  -- parameter denotes the maximal distance two features can be apart without
  -- creating a new Range. Returns an array of Ranges denoting parts the of
  -- <seqid> covered by features.
  function GenomeTools_feature_index:get_coverage(seqid, maxdist)
    assert(seqid)
    local maxdist = maxdist or 0
    local features = self:get_features_for_seqid(seqid)
    local starpos, endpos
    local minstartpos = nil
    local maxendpos = nil
    local ranges = {}
    local coverage = {}

    -- collect all feature ranges
    for i, feature in ipairs(features) do
      ranges[#ranges+1] = feature:get_range()
    end
    -- sort feature ranges
    ranges = gt.ranges_sort(ranges)

    -- compute and store coverage
    for i, range in ipairs(ranges) do
      startpos, endpos = range:get_start(), range:get_end()
      if i == 1 then
        minstartpos = startpos
        maxendpos   = endpos
      else
        -- assert(startpos >= minstartpos)
        if (startpos > maxendpos + maxdist) then
          -- new region started
          coverage[#coverage+1] = gt.range_new(minstartpos, maxendpos)
          minstartpos = startpos
          maxendpos   = endpos
        else
          -- continue old region
          maxendpos = (endpos > maxendpos) and endpos or maxendpos
        end
      end
    end
    -- add last region
    coverage[#coverage+1] = gt.range_new(minstartpos, maxendpos)
    return coverage
  end

  -- Returns an array of Ranges denoting parts of <seqid> which are covered by
  -- at least one marked feature. Internally, get_coverage() is called and the
  -- <maxdist> is passed along.
  function GenomeTools_feature_index:get_marked_regions(seqid, maxdist)
    assert(seqid, "missing seqid argument")
    local coverage = self:get_coverage(seqid, maxdist)
    local marked = {}
    for _,range in ipairs(coverage) do
      local features = feature_index:get_features_for_range(seqid, range)
      if gt.features_contain_marked(features) then
        marked[#marked+1] = range
      end
    end
    return marked
  end

  -- Render to PNG file <png_file> for <seqid> in <range> with optional <width>.
  -- If no <png_file> is given os.tmpname() is called to create one.
  -- Returns name of written PNG file.
  function GenomeTools_feature_index:render_to_png(seqid, range, png_file, width)
    assert(seqid and range)
    png_file = png_file or os.tmpname()
    if not width then width = 1600 end
    local diagram = gt.diagram_new(self, seqid, range)
    local render =  gt.render_new()
    render:to_png(diagram, png_file, width)
    return png_file
  end

  -- Show all sequence IDs.
  function GenomeTools_feature_index:show_seqids()
    for _,seqid in ipairs(feature_index:get_seqids()) do
      print(seqid)
    end
  end

  -- Returns all features from <feature_index>.
  function GenomeTools_feature_index:get_all_features()
    local seqids = self:get_seqids()
    local all_features = {}
    for _, seqid in ipairs(seqids) do
      local seqid_features = self:get_features_for_seqid(seqid)
      for _, feature in ipairs(seqid_features) do
        all_features[#all_features + 1] = feature
      end
    end
    return all_features
  end
end
