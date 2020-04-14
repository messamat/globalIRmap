#### Only run these functions interactively in console, never in RMarkdown document

library(workflowr)

# Tell Git your name and email - Only need to run it once per computer
wflow_start(directory=getwd(), existing=T)

wflow_build()
wflow_publish("analysis/license.Rmd", "First publication")
wflow_status()
