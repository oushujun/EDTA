if not gtdata_doc_dir then gtdata_doc_dir = "./" end
print ([[
File format for option '-rule_files':

The files supplied to option '-rule_files' define a function for
filtering by user given criteria (see example below):
]])
print(io.open(gtdata_doc_dir.."luafilter_function.lua"):read("*a"))
print([[The above function iterates over all children of 'gn' and
checks whether there is a node of type 'exon'. If there is such a
node the function returns 'false', indicating that the parent node
'gn' will not be sorted out.]])
print([[

NOTE:]])
print([[The function must be named 'filter' and must return 'false',
indicating that the node survived the filtering process.]])
