  function filter(gn)
    target = "exon"
    for curnode in gn:children() do
      if (curnode:get_type() == target) then
        return false
      end
    end
    return true
  end
