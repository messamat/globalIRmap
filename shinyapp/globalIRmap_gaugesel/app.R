library(data.table)
library(drake)
library(shiny)
library(leaflet)
library(leafpop)
library(fst)


shinydat <- file.path('shinyapp', 'globalIRmap_gaugesel', 'www') #If run within globalIRmap R project
shinydat <- 'www' #If deployed online
cacheddat <- file.path(shinydat, 'shinygdat.fst')

if (!file.exists(cacheddat)) {
    source('R/IRmapping_packages.R')
    source('R/IRmapping_functions.R')
    
    #Create data.table with gauge id, characteristics, X, Y, reason for removal, and embedded image
    loadd(gaugep)
    loadd(GRDCgaugestats)
    loadd(GSIMgaugestats)
    loadd(path_resdir)
    
    gaugestats_analyzed_format <- analyzemerge_gaugeir(in_GRDCgaugestats = GRDCgaugestats,
                                                       in_GSIMgaugestats = GSIMgaugestats,
                                                       yearthresh = 1800,
                                                       in_gaugep = gaugep,
                                                       inp_resdir = path_resdir,
                                                       plotseries = FALSE)
    
    geocols <- c('X', 'Y')
    
    gaugep_flags <- merge(gaugestats_analyzed_format$flags, gaugep, by=c('GRDC_NO', 'gsim_no'), all.y=T) %>%
        merge(rbind(rbindlist(GRDCgaugestats),
                    rbindlist(GSIMgaugestats), fill=T, use.names=T),
              by=c('GRDC_NO', 'gsim_no'), all.x=T)
    
    gaugep_flags[, (geocols) := as.data.table(st_coordinates(gaugep_flags$geometry))] %>%
        .[, flag := factor(flag, levels=c(NA, 'inspected', 'removed'))] %>%
        .[dor_pc_pva>=500, `:=`(flag = 'removed',
                                 comment = 'excessive degree of regulation')] %>% #Remove those with more than 50% DOR
        .[totalYears_kept_o1800<10, `:=`(flag = 'removed',
                                         comment = '< 10 valid years of data')] %>%
        .[, geometry := NULL] %>%
        .[, GAUGE_NO := fifelse(is.na(gsim_no), GRDC_NO, gsim_no)]
    
    write_fst(gaugep_flags, cacheddat)
} else {
    gaugep_flags <- read_fst(cacheddat) %>% setDT
}
    

pal <- colorFactor(
    palette = "Set1",
    domain = unique(gaugep_flags$flag))

gitplots_rooturl <- "https://github.com/messamat/GSIMires/blob/main/plots/"


gaugep_flags[, `:=`(
    hoverlabs = paste(
        sep = "<br/>",
        fifelse(is.na(gsim_no), paste('GRDC_NO:',GRDC_NO), paste('GSIM_NO:',gsim_no)),
        paste0('Flag: ', fifelse(is.na(flag), 'kept', as.character(flag))),
        paste0('Reason: ',comment)),
    graphurl = paste0(gitplots_rooturl, GAUGE_NO, ".png?raw=true")
)
]

gperkept<- gaugep_flags[(intermittent_o1800 == 0) &
                            (is.na(flag) | flag != 'removed'),]
girkept <- gaugep_flags[(intermittent_o1800 == 1) &
                            (is.na(flag) | flag != 'removed'),]
gperremoved <- gaugep_flags[(intermittent_o1800 == 0 | is.na(intermittent_o1800))
                            & flag == 'removed',]
girremoved <- gaugep_flags[(intermittent_o1800 == 1 | is.na(intermittent_o1800))
                           & flag == 'removed',]

######################### SHINYAPP #############################################
ui <- bootstrapPage(

    tabPanel(
        "World", value="wrld",
        div(class="outer",
            tags$style(
                type = "text/css",
                ".outer {position: fixed; top: 30px; left: 0; right: 0; bottom: 0;
                   overflow: scroll; padding: 20px}", #Padding to make table fully visible
                ".leaflet .legend i{border-radius: 50%; width: 10px;height: 10px;margin-top: 4px;}", #Allow for circle legend
                ".navbar {font-size: 16px}", ".navbar-default .navbar-brand {font-size: 20px}", #Increase font size of navigation bar text
                ".tooltip-inner {text-align: left}"), #Left-align tooltips that display upon hover

            leafletOutput("Worldmap", width = "100%", height = "100%"),
            # absolutePanel(id = "controls", class = "panel panel-default", fixed = F,
            #               draggable = FALSE, top = 60, left = "auto", right = 20, bottom = "auto",
            #               width = 330, height = "auto",
            #               wellPanel(
            #                   p(HTML("<strong>Hydrograph</strong>"))
            #               )
            # ),
            tags$div(id="cite2",
                     HTML('<a href="insert link"
                      target="_blank">Messager et al. (2021) "Global prevalence of non-perenniall rivers and streams - in review"</a>')
            )
        )
    )
)



server <- function(input, output, session) {
    #Custom function to make circle legends from https://stackoverflow.com/questions/37446283/creating-legend-with-circles-leaflet-r
    addLegendCustom <- function(map, colors, labels, sizes, opacity = 0.75){
        colorAdditions <- paste0(colors, "; width:", sizes, "px; height:", sizes, "px")
        labelAdditions <- paste0("<div style='display: inline-block;height: ", sizes, "px;margin-top: 4px;line-height: ", sizes, "px;'>", labels, "</div>")

        return(addLegend(map, colors = colorAdditions, labels = labelAdditions, opacity = opacity, position="bottomleft", title="Streamgages",layerId="gage"))
    }
    #Set up basic map
    output$Worldmap <- renderLeaflet({
        leaflet() %>%
            setView(lng = 34, lat = 28, zoom = 2) %>%
            addTiles(urlTemplate = "//{s}.tiles.mapbox.com/v3/jcheng.map-5ebohr46/{z}/{x}/{y}.png",
                     attribution = 'Maps by <a href="http://www.mapbox.com/">Mapbox</a>',
                     options = tileOptions(opacity=0.5)) %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery",
                             options = providerTileOptions(opacity=0.75)) %>%
            addProviderTiles(providers$Esri.OceanBasemap, group = "World Physical",
                             options = providerTileOptions(opacity=0.75)) %>%
            addCircleMarkers(data= gperkept,
                             weight = 0, radius=3, lng = ~X, lat = ~Y,
                             label = ~lapply(hoverlabs, HTML),
                             popup = popupImage(gperkept$graphurl, src = "remote",
                                                width = 500, height = 500),
                             group='used in modeling - perennial', color='blue',fillOpacity = 0.3) %>%
            addCircleMarkers(data=girkept,
                             weight = 0, radius=3, lng = ~X, lat = ~Y,
                             label = ~lapply(hoverlabs, HTML),
                             popup = popupImage(girkept$graphurl, src = "remote",
                                                width = 500, height = 500),
                             group="used in modeling - non-perennial", color='red',fillOpacity = 0.3) %>%
            addCircleMarkers(data=gperremoved,
                             weight = 0, radius=3, lng = ~X, lat = ~Y,
                             label = ~lapply(hoverlabs, HTML),
                             popup = popupImage(gperremoved$graphurl, src = "remote",
                                                width = 500, height = 500),
                             group="excluded from modeling - perennial or undertermined", color='blue',fillOpacity = 1) %>%
            addCircleMarkers(data=girremoved,
                             weight = 0, radius=3, lng = ~X, lat = ~Y,
                             label = ~lapply(hoverlabs, HTML),
                             popup = popupImage(girremoved$graphurl, src = "remote",
                                                width = 500, height = 500),
                             group="excluded from modeling  - non-perennial", color='red',fillOpacity = 1) %>%
            addLayersControl(
                position="bottomright",
                overlayGroups = c("used in modeling - perennial",
                                  "used in modeling - non-perennial",
                                  "excluded from modeling - perennial or undertermined",
                                  "excluded from modeling  - non-perennial"),
                options = layersControlOptions(collapsed = TRUE)) %>%
            addLegend(colors=c('blue', 'red'),
                      position="topright",
                      labels=c("excluded from modeling  - perennial",
                               "excluded from modeling  - non-perennial"),
                      opacity = 1
            ) %>%
            addLegend(colors=c('blue', 'red'),
                      position="topright",
                      labels=c("used in modeling  - perennial",
                               "used in modeling - non-perennial"),
                      opacity = 0.3
            )
    })
}


shinyApp(ui, server)
