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

-- Returns true if the given array of <features> contains a marked feature,
-- false otherwise.
function features_contain_marked(features)
  assert(features)
  for _, feature in ipairs(features) do
    if feature:contains_marked() then
      return true
    end
  end
  return false
end

-- Print the given array of <features> to stdout.
function features_show(features)
  assert(features)
  local gff3_visitor = gt.gff3_visitor_new()
  for _, features in ipairs(features) do
    features:show(gff3_visitor)
  end
end

-- Return all marked <features> (an array) as an array or nil if <features>
-- contains no marked features.
function features_get_marked(features)
  assert(features)
  local marked_features = nil
  if features_contain_marked(features) then
    marked_features = {}
    for _, feature in ipairs(features) do
      if feature:contains_marked() then
        local gni = gt.genome_node_iterator_new(feature)
        local node = gni:next()
        while node do
          if node:is_marked() then
            marked_features[#marked_features + 1] = node
          end
          node = gni:next()
        end
      end
    end
  end
  return marked_features
end

-- Print all marked <features> (an array) to stdout.
function features_show_marked(features)
  assert(features)
  if features_contain_marked(features) then
    for _, feature in ipairs(features) do
      feature:show_marked()
    end
  end
end

local function create_gene_from_mRNA(mRNA)
  assert(mRNA)
  assert(mRNA:get_type() == "mRNA")
  local gene = gt.genome_feature_new(mRNA:get_seqid(), "gene", mRNA:get_range(),
                                     mRNA:get_strand())
  gene:set_source(mRNA:get_source())
  local gni = gt.genome_node_iterator_new_direct(mRNA)
  local old_child = gni:next()
  while (old_child) do
    local new_child = gt.genome_feature_new(old_child:get_seqid(),
                                            old_child:get_type(),
                                            old_child:get_range(),
                                            old_child:get_strand())
    new_child:set_source(old_child:get_source())
    gene:is_part_of_genome_node(new_child)
    old_child = gni:next()
  end
  return gene
end

-- Return an array of genome features which contains a separate gene feature for
-- each mRNA in <in_features>.
function features_mRNAs2genes(in_features)
  assert(in_features)
  local out_features = {}
  for _, in_feature in ipairs(in_features) do
    if in_feature:get_type() == "gene" then
      local gni = gt.genome_node_iterator_new_direct(in_feature)
      local child = gni:next()
      while (child) do
        if child:get_type() == "mRNA" then
          out_features[#out_features + 1] = create_gene_from_mRNA(child)
        end
        child = gni:next()
      end
    end
  end
  return out_features
end

-- Return an array with the sequences of the given features.
function features_extract_sequences(features, type, join, region_mapping)
  local sequences = {}
  for _, feature in ipairs(features) do
    local sequence = feature:extract_sequence(type, join, region_mapping)
    if type then
      sequences[#sequences + 1] = sequence
    end
  end
  return sequences
end
