#Code based on https://github.com/matthewstrasiotto/mandrake/blob/master/vignettes/mandrake.Rmd

source('_drake.R')

drake::vis_drake_graph(plan, targets_only = T)


# Load the column specification for mandrake.
cache <- drake::get_cache()
lookup_cache <- mandrake::load_package_colspec("mandrake")

# Attach the column documentation to our plan.
plan_extracted <- plan %>%
  mandrake::decorate_plan(cache,
                          #group = "cluster_id",
                          lookup_cache = lookup_cache)

graph_info <- drake_graph_info(
  plan_extracted,
  hover = T,
  # group = "cluster_id",
  # clusters = c("summ", "coef"),
  build_times = "none",
  on_select_col = "desc"
)

graph <- render_drake_graph(
  graph_info,
  on_select = "embedHandler",
  ncol_legend = 4
) %>% mandrake::attach_dependencies()

graph %>%
  visNetwork::visHierarchicalLayout(
    direction = "LR"
  ) %>%
  visNetwork::visEdges(
    smooth = list(type = "cubicBezier", forceDirection = "horizontal")
  )
