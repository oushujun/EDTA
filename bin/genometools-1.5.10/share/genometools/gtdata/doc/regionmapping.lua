if not gtdata_doc_dir then gtdata_doc_dir = "./" end
print ([[

File format for option '-regionmapping':

The file supplied to option -regionmapping defines a ``mapping''.  A mapping
maps the `sequence-region` entries given in the 'GFF3_file' to a sequence file
containing the corresponding sequence. Mappings can be defined in one of the
following two forms:
]])
print(io.open(gtdata_doc_dir.."regionmapping_table.lua"):read("*a"))
print([[
or
]])
print(io.open(gtdata_doc_dir.."regionmapping_function.lua"):read("*a"))
print([[
The first form defines a Lua (http://www.lua.org) table named ``mapping''
which maps each sequence region to the corresponding sequence file.
The second one defines a Lua function ``mapping'', which has to return the
sequence file name when it is called with the `sequence_region` as argument.]])
