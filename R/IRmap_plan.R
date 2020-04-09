plan <- drake_plan(
  filestructure = def_filestructure(),
  gaugep = read_gaugep(in_filestructure=filestructure, dist=200),
  gauged_filenames = read_gauged_paths(filestructure, gaugep),
  gaugestats = future_map(gauged_filenames, #To map
                          comp_durfreq, #Function to run on each file name
                          maxgap=20, monthsel=NULL, #Other arguments in function
                          .progress = TRUE),
  gaugestats_format = format_gaugestats(gaugestats, gaugep),
  predvars = selectformat_predvars(filestructure, gaugestats_format)
)

  