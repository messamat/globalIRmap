#### Only run these functions interactively in console, never in RMarkdown document

library(workflowr)

# Tell Git your name and email - Only need to run it once per computer
wflow_start(directory=getwd(), existing=T)

wflow_build()
wflow_publish('analysis/index.Rmd')
wflow_status()

wflow_publish(c("analysis/_site.yml", "analysis/diagnostics.Rmd"),
               "Add tab to website, adjust figure sizes")

#To add new tab, add
