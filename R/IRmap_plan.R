plan <- drake_plan(
  file_structure = def_filestructure(),
  gaugep = read_gaugep(file_structure, dist=200),
  gauged_filenames = read_gauged_paths(file_structure, gaugep),
  gaugestats = target(durfreq_indiv(gauged_filenames, maxgap=20, monthsel=NULL), 
                       dynamic = map(gauged_filenames))
)
  