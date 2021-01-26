library(data.table)
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


    analyzemerge_gaugeir <- function(in_GRDCgaugestats, in_GSIMgaugestats, yearthresh,
                                     in_gaugep, inp_resdir, plotseries = FALSE) {
        ### Analyze GSIM data ########################################################
        GSIMstatsdt <- rbindlist(in_GSIMgaugestats)

        #------ Remove stations with unstable intermittent flow regime
        #Remove those which have at least one day per year of zero-flow day but instances
        #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
        GSIMtoremove_unstableIR <- data.table(
            gsim_no = GSIMstatsdt[(mDur_o1800 >= 1) & (!movinginter_o1800), gsim_no],
            flag = 'removed',
            comment = 'automatic filtering: at least one no-flow day/year on average but no zero-flow event during >= 20 years'
        )

        #------ Remove stations based on examination of plots and data series
        #Outliers from examining plots of ir time series (those that were commented out were initially considered)
        GSIMtoremove_irartifacts <- list(
            c('AR_0000014', 'removed', "large gaps in data, changed flow permanence"),
            c('AT_0000021', 'removed', "single flow intermittency event, probably gap in data"),
            c('AT_0000026', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
            c('AT_0000038', 'removed', "large gaps in data, single flow intermittency event at the end"),
            c('AT_0000059', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
            c('BR_0000286', 'removed', "Tapajos river. Impossible that it dries out."),
            c('BR_0000557', 'inspected', "Confirmed dry channel visually on satellite imagery"),
            c('BR_0000581', 'removed', "Single flow intermittency event at the end"),
            c('BR_0000662', 'removed', "flow regulated. unsure when it started, large data gap"),
            c('BR_0000664', 'inspected', "change of regime in last few years due to regulation, but originally IRES"),
            c('BR_0000706', 'removed', "too many data gaps to determine long-term flow permanence"),
            c('BR_0000717', 'removed', "flow regulated, changed from perennial to IRES"),
            c('BR_0000778', 'removed', "too many data gaps to determine long-term flow permanence"),
            c('BR_0000786', 'removed', "only one flow intermittency event"),
            c('BR_0000862', 'removed', "only one flow intermittency event, too many data gaps"),
            c('BR_0001011', 'removed', "only one flow intermittency event, too many data gaps"),
            c('BR_0001104', 'removed', "only one flow intermittency event at the end"),
            c('BR_0001116', 'removed', "seems to have changed flow permanence, many gata gaps"),
            c('BR_0001133', 'removed', "changed flow permanence"),
            c('CA_0001057', 'removed', "only one flow intermittency event in 28 years"),
            c('CA_0003473', 'inspected', "only one flow intermittency event but outlet of natural lake"),
            c('CA_0003488', 'removed', "only one flow intermittency event"),
            c('CA_0003526', 'removed', "only one flow intermittency event over the winter"),
            c('CA_0003544', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000004', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000009', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000010', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000012', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000013', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000004', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000021', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000022', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000026', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000029', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000038', 'removed', "only one flow intermittency event, temporary dewatering by dam"),
            c('CN_0000043', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000047', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('CN_0000062', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('ES_0000078', 'removed', "changed flow permanence, first 20 years with 0 flow"),
            c('ES_0000087', 'removed', "0 flow seems driven by exceptional drought, not throughout record"),
            c('ES_0000388', 'removed', "too many data gaps to determine long-term flow permanence; and looks regulated"),
            c('ES_0000444', 'removed', "only one flow intermittency event"),
            c('ES_0000525', 'removed', "very discontinued record. Seems to have changed flow permanence"),
            c('ES_0000581', 'removed', "appears unreliable"),
            c('ES_0000676', 'inspected', "flow regulated, but originally IRES"),
            c('ES_0000729', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('ES_0000784', 'removed', "same as ES_0000785"),
            c('ES_0000785', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('ES_0000786', 'removed', "abrupt decreases to 0 flow, probably gaps in data"),
            c('ES_0000794', 'removed', "same as ES_0000832"),
            c('ES_0000816', 'removed', "abrupt decreases to 0 flow"),
            c('ES_0000818', 'removed', "same as ES_0000785 and ES_0000784, which are not intermittent"),
            c('ES_0000856', 'removed', "only one flow intermittency event"),
            c('ES_0000892', 'removed', "appears perennial with only one period of flow intermittency, data gaps"),
            c('ES_0000906', 'removed', "abrupt decreases to 0 flow"),
            c('ES_0000910', 'removed', "only one flow intermittency event"),
            c('ES_0000910', 'removed', "only one flow intermittency event"),
            c('ES_0000958', 'removed', "only one flow intermittency event, abrupt decrease"),
            c('ES_0000986', 'removed', "only one flow intermittency event, abrupt decrease"),
            c('ES_0000996', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('ES_0001020', 'removed', "too many data gaps to determine long-term flow permanence"),
            c('ES_0001052', 'removed', "regulated, changed flow permanence"),
            c('ES_0001082', 'removed', "regulated, changed flow permanence"),
            c('ES_0001085', 'removed', "only one flow intermittency event at the end"),
            c('ES_0001116', 'removed', "only one year of flow intermittency events at the end"),
            c('ES_0001162', 'removed', "abrupt decreases to 0 flow"),
            c('FI_0000107', 'removed', "lake inlet, maybe just become standing water when high water level"),
            c('FI_0000156', 'removed', "no 0 flow for first 60% of record"),
            c('FR_0000052', 'removed', "no 0 flow for first 20 years of record"),
            c('HU_0000017', 'removed', "only one 0 flow occurence in 25 years"),
            c('IN_0000014', 'removed', "no 0 flow occurrence in first 25 years"),
            c('IN_0000023', 'removed', "no 0 flow occurrence in first 20 years"),
            c('IN_0000045', 'removed', "no 0 flow occurrence in first 20 years"),
            c('IN_0000046', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000050', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000063', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000064', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000105', 'removed', "regulated -- no records pre-1968 time of dam building"),
            c('IN_0000113', 'removed', "only one 0 flow occurence"),
            c('IN_0000124', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000125', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000134', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000159', 'removed', "regulated -- no records pre-1968 time of dam building"),
            c('IN_0000170', 'removed', "only one 0 flow occurence"),
            c('IN_0000190', 'inspected', "looked regulated but isn't in the end"),
            c('IN_0000255', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000105', 'removed', "regulated -- no records pre-1968 time of dam building"),
            c('IN_0000280', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('IN_0000309', 'inspected', "before building of reservoir in 1988"),
            c('IN_0000312', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('MZ_0000010', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('NA_0000050', 'removed', "unreliable record, data gaps and interpolated values"),
            c('NO_0000018', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('NO_0000028', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('NO_0000030', 'removed', "only one 0 flow occurence"),
            c('NO_0000044', 'removed', "0 flow values only at the beginning, probably data gaps"),
            c('NO_0000090', 'removed', "0 flow values at the beginning probably data gaps, otherwise only one 0 flow event"),
            c('RU_0000189', 'removed', "only one 0 flow occurence"),
            c('RU_0000265', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('RU_0000358', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('SE_0000058', 'removed', "downstream of dam, changed flow permanence. and experiences ice"),
            c('TZ_0000018', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('US_0005106', 'removed', "probably changed flow permanence from perennial to non-perennial"),
            c('US_0005874', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('US_0005876', 'removed', "probably changed flow permanence from perennial to non-perennial"),
            #'US_0001855', #maybe rounded values --- checked
            #'US_0001861', #maybe rounded values --- checked
            #'US_0001868', #maybe rounded avlues --- checked
            c('US_0002247', 'removed', "regulated, GRDC 4149415 upstream not intermittent"),
            c('US_0002248', 'removed', "regulated, just downstream of US_0002247"),
            c('US_0002791', 'removed', "on usgs website: close proximity to the Ohio River. During periods of high water on the Ohio River, the computed discharge at this site may be incorrectly displayed due to the backwater effect created."),
            #'US_0003591', #maybe rounded values --- checked
            # 'US_0003774', #maybe rounded values --- checked
            # 'US_0003836', #maybe rounded values --- checked
            # 'US_0004023', #maybe rounded values --- checked
            # 'US_0004216', #maybe rounded values --- checked
            # 'US_0004232', #maybe rounded values --- checked
            # 'US_0004658', #maybe rounded values --- checked
            c('US_0004773', 'removed', "regulated, no records prior to reservoir buiding"),
            c('US_005099', 'inspected', "confirmed dry bed on satellite imagery"),
            # 'US_0005161', #maybe rounded values --- checked
            # 'US_0005177', #maybe rounded values --- checked
            c('US_0005303', 'inspected', "looks fine on usgs website and imagery. just small"),
            c('US_0005596', 'inspected', "looks fine on usgs website and imagery. just small"),
            c('US_0005597', 'inspected', "looks fine on usgs website and imagery. just small"),
            # 'US_0005622', #maybe rounded values --- checked
            # 'US_0005623', #maybe rounded values --- checked
            # 'US_0005684', #maybe rounded values --- checked
            # 'US_0005687', #maybe rounded values --- checked
            c('US_0005732', 'removed', "regulated by Lake Arcadia reservoir, changed flow permanence"),
            #'US_0005859', #maybe rounded values --- checked
            #'US_0005879', #maybe rounded values --- checked
            # 'US_0006073', #maybe rounded values --- checked
            c('US_0006103', 'removed', "regulated, changed flow permanence"),
            # 'US_0006109', #maybe rounded values --- check
            # 'US_0006154', #maybe rounded values --- check
            c('US_0006155', 'removed', "regulated, did not change flow permanence but same as 0006156"),
            c('US_0006156', 'inspected', "regulated, but did not change flow permanence"),
            c('US_0006206', 'inspected', "looks fine on usgs website and imagery"),
            #'US_0006301', #maybe rounded values --- checked
            #'US_0006327', #maybe rounded values --- checked
            #'US_0006387', #maybe rounded values --- checked
            c('US_0006396', 'removed', "erroneous values pre-1960"),
            c('US_0006440', 'removed', "probably changed flow permanence from perennial to non-perennial"),
            c('US_0006537', 'removed', "changed flow permanence from perennial to non-perennial"),
            c('US_0006574', 'inspected', "suspected regulation, but doesn't seem to be the case"),
            #'US_0006975', #maybe rounded values --- checked
            #'US_0006984', #maybe rounded values --- checked
            #'US_0006985', #maybe rounded values --- checked
            #'US_0006986', #maybe rounded values --- checked
            c('US_0008607', 'inspected', "regulated since 1978, schanged flow permanence from perennial to non-perennial"),
            c('US_0008687', 'inspected', "regulated since 1935, schanged flow permanence from perennial to non-perennial"),
            # 'US_0008726', #maybe rounded values --- checked
            # 'US_0008779', #maybe rounded values --- checked
            c('ZA_0000074', 'removed', "abrupt decrease to 0 flow, probably gaps in data"),
            c('ZA_0000268', 'removed', "probably changed flow permanence"),
            c('ZA_0000084', 'removed', "probably changed flow permanence"),
            c('ZA_0000270', 'removed', "changed flow permanence")
        ) %>%
            do.call(rbind, .) %>%
            as.data.table %>%
            setnames(c('gsim_no', 'flag', 'comment'))

        # readformatGSIMmon(
        #   GSIMstatsdt[gsim_no == 'US_0000546', path]) %>%
        #   .[MIN != 0, min(MIN)]
        # sort(GSIMtoremove_unstableIR)
        # 'ES_0000806' %in% GSIMstatsdt$gsim_no
        # GSIMstatsdt[gsim_no == 'ES_0000818', ]
        #
        #-----  Check flags in winter IR for GSIM
        wintergaugesall_GSIM <- plot_winterir(
            dt = GSIMstatsdt, dbname = 'gsim', inp_resdir = inp_resdir,
            yearthresh = 1800, plotseries = plotseries)
        #
        # #Check suspicious canadian ones
        # GSIMwintermeta <- in_gaugep[in_gaugep$gsim_no %in% wintergaugesall_GSIM$gsim_no,]
        # canadians_toinspect <- in_gaugep[in_gaugep$gsim_no %in%
        #                                    paste0('CA_000', c(3469, 3473, 3526, 3544, 6082, 6122)),]$reference_no
        #
        # if (!dir.exists(hy_dir())) download_hydat()
        # cancheck <- lapply(canadians_toinspect, function(refno) {
        #   merge(hy_daily(station_number = refno),
        #         hy_stn_regulation(station_number = refno),
        #         by='STATION_NUMBER') %>%
        #     setDT
        # }) %>%
        #   rbindlist
        # cancheck[REGULATED==T, unique(STATION_NUMBER)] #No regulated station
        # cancheck[Value==0, .N, by=.(STATION_NUMBER, Symbol)]

        #E - Estimate:  no measured data available for the day or missing period,
        #    and the water level or streamflow value was estimated by an indirect method
        #A - Partial Day:  daily mean value was estimated despite gaps of more than
        #    120 minutes in the data string or missing data not significant enough to
        #    warrant the use of the E symbol.
        #B - Ice conditions: value was estimated with consideration for the presence
        #    of ice in the stream. Ice conditions alter the open water relationship
        #    between water levels and streamflow.
        #D - Dry: stream or lake is "dry" or that there is no water at the gauge.
        #    This symbol is used for water level data only.
        #R - Revised: The symbol R indicates that a revision, correction or addition
        # `  has been made to the historical discharge database after January 1, 1989.

        # ggplot(cancheck[Value > 0, ], aes(x=Date, y=Value, color=Symbol)) +
        #   geom_vline(data=cancheck[is.na(Value),], aes(xintercept = Date), color='grey', alpha=1/4) +
        #   geom_point(alpha=1/6) +
        #   geom_point(data=cancheck[Value==0,]) +
        #   facet_wrap(~STATION_NUMBER, scales='free') +
        #   theme_classic()
        #
        #Check others 'CN_0000047', 'NO_0000018', 'RU_0000089',
        #'RU_0000391', 'RU_0000393', 'RU_00000395', 'RU_0000436',
        #'RU_0000470', 'US_0008687')
        #check <- readformatGSIMmon(GSIMstatsdt[gsim_no == 'US_0008687',path])

        #Remove
        GSIMtoremove_winterIR <- list(
            c('CA_0003473', 'removed', "sudden peak — unsure about estimated discharge under ice conditions"),
            c('CN_0000047', 'removed', "Anomalous change from near 0 discharge to 150 m3/s, no explanation"),
            c('RU_0000391', 'removed', "Stopped recording during the winter the last ~10 years. Maybe questionable winter data"),
            c('RU_0000393', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data"),
            c('RU_0000436', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data"),
            c('RU_0000470', 'removed', "Didn't record during the winter for the first 20 years. Maybe questionable winter data")
        ) %>%
            do.call(rbind, .) %>%
            as.data.table %>%
            setnames(c('gsim_no', 'flag', 'comment'))

        #-----  Check flags in coastal IR for GSIM
        GSImcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GSIMstatsdt,
                                           dbname = 'gsim', inp_resdir = inp_resdir,
                                           yearthresh = 1800, plotseries = plotseries)
        #Already removed suspect ones

        ### Analyze GRDC data ########################################################################
        GRDCstatsdt <- rbindlist(in_GRDCgaugestats)

        #Remove all gauges with 0 values that have at least 99% of integer values as not reliable (see GRDC_NO 6140700 as example)
        GRDCtoremove_allinteger <- data.table(
            GRDC_NO = GRDCstatsdt[integerperc_o1800 >= 0.99 &
                                      intermittent_o1800 == 1, GRDC_NO],
            flag = 'removed',
            comment = 'All integer discharge values'
        )

        #Remove those which have at least one day per year of zero-flow day but instances
        #of no zero-flow day within a 20-year window — except for three gauges that have a slight shift in values but are really IRES
        GRDCtoremove_unstableIR <- data.table(
            GRDC_NO = GRDCstatsdt[(mDur_o1800 >= 1) & (!movinginter_o1800) &
                                      !(GRDC_NO %in% c(1160115, 1160245, 4146400)), GRDC_NO],
            flag = 'removed',
            comment = 'automatic filtering: at least one no-flow day/year on average but no zero-flow event during >= 20 years'
        )

        #Outliers from examining plots of ir time series (those that were commented out were initially considered)
        GRDCtoremove_irartifacts <- list(
            c(1134300, 'removed', 'changed flow permanence from perennial to non-perennial, large data gaps'),
            c(1134500, 'removed', 'only 1 occurrence of 0 flow values'),
            c(1159302, 'removed', 'abrupt decreases to 0 flow values'),
            c(1159302, 'removed', 'unreliable record, isolated 0s, sudden jumps and capped at 77'),
            c(1159320, 'inspected', "0 values for the first 14 years but still apparently originally IRES"),
            c(1159325, 'inspected', "0 values for most record. probably episodic and due to series of agricultural ponds"),
            c(1159510, 'inspected', "0 values for most record, on same segment as 1159511 but seems unreliable"),
            c(1159520, 'inspected', "values seem capped after 1968, otherwise seem fine. Could just be rating curve"),
            c(1159830, 'removed', 'only one occurence of 0 flow values'),
            c(1160101, 'removed', 'abrupt decreases to 0 flow values'),
            c(1160340, 'removed', 'abrupt decreases to 0 flow values'),
            c(1160378, 'inspected', 'appears IRES before regulation, reservoir fully dry on satellite imagery, so keep as intermittent even when not regulated'),
            c(1160420, 'removed', 'decrease to 0 appears a bit abrupt but ok'),
            c(1160435, 'removed', 'unreliable record, abrupt decreases to 0, capped'),
            c(1160470, 'removed', 'unreliable record, probably change of rating curve in 1947, mostly missing data until 1980 but truly intermittent based on imagery'),
            c(1160540, 'inspected', 'only 0 - nodata for first 15 years. seemingly good data post 1979 and IRES'),
            c(1160635, 'inspected', '0 values post-1985 seem abrupt but enough values otherwise to make it intermittent'),
            c(1160670, 'removed', 'regulated'),
            c(1160675, 'removed', 'some outliers but otherwise most 0 values seem believable'),
            c(1160780, 'removed', 'unreliable record, abrupt decreases to 0 flow values, large data gaps'),
            c(1160785, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1160795, 'removed', 'abrupt decreases to 0 flow values'),
            c(1160800, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1160840, 'removed', 'only 2 zero flow values are believable, others are outliers'),
            c(1160850, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1160880, 'removed', 'unreliable record. Tugela river, perennial'),
            c(1160881, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1160900, 'removed', 'most 0 values look like outliers, abrupt decreases'),
            c(1160911, 'removed', 'most 0 values look like outliers, abrupt decreases'),
            c(1160971, 'removed', 'most 0 values look like outliers, abrupt decreases'),
            c(1160975, 'removed', 'most 0 values look like outliers, abrupt decreases'),
            c(1196102, 'removed', 'unreliable record, large data gaps, hard to tell original flow permanence'),
            c(1196141, 'removed', "doesn't look reliable, hard to assess long term flow permanence"),
            c(1196160, 'inspected', 'some outlying 0 flow values but most are good'),
            c(1197500, 'removed', 'only one flow intermittency event, abrupt decrease to 0'),
            c(1197540, 'removed', 'abrupt decreases to 0'),
            c(1197591, 'removed', 'abrupt decreases to 0'),
            c(1197700, 'removed', 'abrupt decreases to 0'),
            c(1197740, 'removed', 'some outlying 0 flow values but most are good'),
            c(1199100, 'removed', 'most 0 values look like outliers'),
            c(1199200, 'removed', 'abrupt decreases to 0'),
            c(1199410, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1259500, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1259800, 'removed', '0 values come from integer-based part of the record'),
            c(1286690, 'removed', 'changed flow permanence, record too short to determine original flow permanence'),
            c(1259800, 'removed', 'changed flow permanence, only one 0 flow value post 1963'),
            c(1428400, 'removed', '0 values come from integer-based part of the record'),
            c(1428500, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1434200, 'removed', 'almost all integers'),
            c(1434300, 'removed', 'almost all integers'),
            c(1434810, 'removed', '0 values come from integer-based part of the record'),
            c(1491815, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(1491870, 'removed', 'some outlying 0 flow values but most are good'),
            c(1494100, 'removed', 'abrupt decreases to 0'),
            c(1494100, 'inspected', 'regulated but naturally intermittent'),
            c(1591110, 'removed', "doesn't look reliable, changed flow permanence from perennial to non-perennial"),
            c(1591730, 'removed', 'abrupt decreases to 0'),
            c(1733600, 'removed', '0 values come from integer-based part of the record and outliers'),
            c(1837410, 'removed', 'abrupt decreases to 0'),
            c(1837430, 'inspected', 'nearly same as 1837410. Naturally intermittent before dam'),
            c(1897550, 'removed', 'abrupt decreases to 0'),
            c(1898501, 'removed', 'abrupt decreases to 0'),
            c(1992400, 'removed', 'most 0 values look like outliers'),
            c(2588500, 'removed', 'abrupt decreases to 0'),
            c(2588551, 'removed', 'abrupt decreases to 0'),
            c(2588630, 'removed', 'abrupt decreases to 0'),
            c(2588640, 'removed', 'abrupt decreases to 0'),
            c(2588708, 'removed', 'abrupt decreases to 0'),
            c(2588820, 'removed', 'abrupt decreases to 0'),
            c(2589230, 'removed', 'abrupt decreases to 0'),
            c(2589370, 'removed', 'abrupt decreases to 0'),
            c(2591801, 'removed', 'abrupt decreases to 0'),
            c(2694450, 'removed', 'abrupt decreases to 0'),
            c(2969081, 'removed', 'abrupt decreases to 0'),
            c(2999920, 'removed', '0 values come from integer-based part of the record'),
            c(3650460, 'inspected', 'some outlying 0 flow values but most are good'),
            c(3650470, 'removed', '0 values come from integer-based part of the record'),
            c(3650610, 'removed', 'integers pre-1960s but still intermittent after'),
            c(3650640, 'removed', 'abrupt decreases to 0'),
            c(3650649, 'inspected', 'change of flow regime due to dam building but intermittent before'),
            c(3650690, 'inspected', 'abrupt decreases to 0'),
            c(3650928, 'removed', 'most 0 values look like outliers'),
            c(3652050, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(3652135, 'removed', 'only one valid 0-flow event'),
            c(3652200, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(3844460, 'removed', 'abrupt decreases to 0'),
            c(3844460, 'removed', 'abrupt decreases to 0'),
            c(4101451, 'inspected', 'station downstream also has 0s'),
            c(4103700, 'removed', '0 values come from integer-based part of the record'),
            c(4150605, 'inspected', 'just downstream of lwesville dam in Dallas. previously intermittent as well but will be removed anyways as >50% dor'),
            c(4151513, 'inspected', 'looks regulated but will be removed as > 50% regulated'),
            c(4208195, 'removed', 'unreliable record, 0 flow values stem from interpolation'),
            c(4208372, 'removed', 'abrupt decreases to 0 flow, probably data gaps'),
            c(4208655, 'removed', 'insufficient data to tell flow permanence'),
            c(4208855, 'removed', 'insufficient data to tell flow permanence'),
            c(4208857, 'removed', 'only 1 occurrence of 0 flow values'),
            c(4213566, 'removed', 'only 1 occurrence of 0 flow values'),
            c(4213905, 'removed', "regulated, changed flow permanence"),
            c(4214075, 'removed', '0 flow values are data gaps'),
            c(4214297, 'removed', 'only 1 occurrence of 0 flow values'),
            c(4214298, 'removed', 'only 1 occurrence of 0 flow values'),
            c(4234300, 'removed', "regulated, changed flow permanence"),
            c(4243610, 'removed', "regulated, abrupt decrease to 0 probably due to reservoir filling/construction"),
            c(4351710, 'removed', '0 values come from integer-based part of the record and outliers'),
            c(4355500, 'removed', "regulated, outlier 0 flow values"),
            c(4357510, 'removed', "single flow intermittency event, probably gap in data"),
            c(4769200, 'removed', 'only 1 occurrence of 0 flow values'),
            c(4773050, 'removed', 'abrupt decreases to 0'),
            c(5101020, 'removed', "single flow intermittency event, probably gap in data"),
            c(5101101, 'removed', "single flow intermittency event, probably gap in data"),
            c(5101130, 'removed', 'abrupt decreases to 0'),
            c(5101201, 'removed', 'unreliable record, large data gap as 0 flow values'),
            c(5101290, 'removed', '0 flow values before 2000 are outliers, changed flow permanence'),
            c(5101305, 'removed', 'most 0 values look like outliers'),
            c(5101380, 'removed', 'abrupt decreases to 0'),
            c(5109200, 'removed', 'unreliable record, interpolation, large data gap as 0 flow values'),
            c(5109230, 'removed', 'abrupt decreases to 0'),
            c(5202140, 'removed', 'abrupt decreases to 0'),
            c(5202145, 'removed', 'most 0 values look like outliers'),
            c(5202228, 'removed', 'maybe regulated, unreliable record post 1983 accounts for 0 flow values'),
            c(5204170, 'removed', 'changed flow permanence'),
            c(5302251, 'removed', 'large data gap as 0 flow values, otherwise only one flow intermittency event'),
            c(5302261, 'removed', 'large data gap as 0 flow values'),
            c(5405095, 'removed', 'changed flow permanence from perennial to non-perennial'),
            c(5608100, 'removed', 'large data gap as 0 flow values'),
            c(5708200, 'removed', 'changed flow permanence'),
            c(5803160, 'removed', 'large data gap as 0 flow values'),
            c(5864500, 'removed', 'rounded to 10L/s'),
            c(5870100, 'removed', 'rounded to 100L/s'),
            c(6119100, 'removed', 'rounded to 10L/s'),
            c(6442300, 'removed', '0 values come from integer-based part of the record and outliers'),## perfect example of what an integer-based record involves
            c(6444250, 'removed', '0 values come from integer-based part of the record and outliers'),
            c(6444350, 'removed', '0 values come from integer-based part of the record and outliers'),
            c(6444400, 'removed', 'abrupt decreases to 0'),
            c(6935570, 'removed', 'rounded to 10L/s')
        ) %>%
            do.call(rbind, .) %>%
            as.data.table %>%
            setnames(c('GRDC_NO', 'flag', 'comment'))

        #### Check intermittent record
        # checkno <- 6444400 #GRDC_NO
        # check <- checkGRDCzeroes( #Check area around 0 values
        #   GRDCstatsdt, in_GRDC_NO=checkno, period=15, yearthresh=1800,
        #   maxgap=20, in_scales='free', labelvals = F)
        # checkno %in% GRDCtoremove_allinteger #Check whether all integers
        # in_gaugep[in_gaugep$GRDC_NO==checkno & !is.na(in_gaugep$GRDC_NO), "dor_pc_pva"] #check DOR
        # GRDCstatsdt[GRDC_NO == checkno, integerperc_o1800] #Check % integers

        #Outliers from examining plots of perennial time series (those that were commented out were initially considered)
        #Try to find those:
        # whose low flow plateaus could be 0s
        # whose perennial character is dam-driven or maybe irrigation driven (changed from IR to perennial but hard to find)
        # whose missing data are actually 0s
        # whose quality is too low to be reliable
        GRDCtoremove_pereartifacts <- list(
            c(1159800, 'removed', 'regulated'),
            c(1160331, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
            c(1160788, 'removed', 'low-flow plateaus may be overestimated 0 values'),
            c(1593100, 'inspected', 'bad quality but clearly not IRES'),
            c(1593751, 'inspected', 'missing values may contain intermittency, strange regular patterns, maybe interpolated'),
            c(3628200, 'removed', 'appears to change flow permanence'),
            c(3652030, 'removed', '0s in missing years and low flows in other years may also be 0s'),
            c(4101200, 'removed', 'low-flow plateaus are likely overestimated 0 values'),
            c(4115225, 'removed', 'regulated, may have been IRES otherwise'),
            c(4118850, 'removed', 'low-flow plateaus may be overestimated 0 values'),
            c(4125903, 'removed', 'regulated, may have been IRES otherwise'),
            c(4126351, 'removed', 'regulated, may have been IRES otherwise'),
            c(4148850, 'removed', 'many 0 flow values in missing years'),
            c(4151801, 'removed', 'regulated, may have been IRES otherwise, Rio Grande'),
            c(4152651, 'removed', 'regulated by blue mesa reservoir, may have been IRES otherwise, Rio Grande'),
            c(4208610, 'removed', 'too much missing data but if not would be IRES'),
            c(4213055, 'removed', 'too much missing data but if not would be IRES'),
            c(4213802, 'removed', 'identical to 4213801'),
            c(4214320, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
            c(4362100, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
            c(5606090, 'removed', 'low-flow plateaus may be overestimated 0 values'),
            c(5606414, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
            c(6123630, 'removed', 'low-flow plateaus may be overestimated 0 values, missing years have 0 flows'),
            c(6335020, 'removed', 'identical to 6335060'),
            c(6335050, 'removed', 'identical to 6335060'),
            c(6337503, 'removed', 'regulated, cannot tell whether may have been intermittent before'),
            c(6442100, 'removed', 'identical to 6442600'),
            c(6935146, 'removed', 'identical to 6935145'),
            c(6335050, 'removed', 'identical to 6335060'),
            c(6935600, 'removed', 'identical to 6935145')
        ) %>%
            do.call(rbind, .) %>%
            as.data.table %>%
            setnames(c('GRDC_NO', 'flag', 'comment'))

        #---------- Check flags in winter IR
        plot_winterir(dt = GRDCstatsdt, dbname = 'grdc', inp_resdir = inp_resdir,
                      yearthresh = 1800, plotseries = plotseries)
        #Checked for seemingly anomalous 0s. Sudden decreases.
        #Check for flags, check satellite imagery, station name, check for construction of reservoir

        #------ Check time series of stations within 3 km of seawater
        GRDCcoastalirall <- plot_coastalir(in_gaugep = in_gaugep, dt = GRDCstatsdt,
                                           dbname = 'grdc', inp_resdir = inp_resdir,
                                           yearthresh = 1800, plotseries = plotseries)
        #GRDCcoastalirall[, unique(readformatGRDC(path)$Flag), by=GRDC_NO]
        #Nothing obviously suspect beyond those that ad already been flagged

        #Inspect statistics for 4208857, 4213531 as no flow days occurred only one year
        # ID = '6976300'
        # GRDCstatsdt[GRDC_NO == ID,]
        # check <- readformatGRDC(GRDCstatsdt[GRDC_NO == ID,path])
        # unique(check$Flag)
        #
        # plotGRDCtimeseries(GRDCstatsdt[GRDC_NO == ID,], outpath=NULL)


        ### Summarize removal ########################################################################

        #Before cleaning
        GRDCflags <- rbindlist(list(GRDCtoremove_allinteger,
                                    GRDCtoremove_unstableIR,
                                    GRDCtoremove_irartifacts,
                                    GRDCtoremove_pereartifacts
        ))

        GRDCtoremove_all <- GRDCflags[flag=='removed', GRDC_NO]

        GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10, .N]
        GRDCstatsdt[intermittent_o1800 == 1 & totalYears_kept_o1800 >= 10 &
                        !(GRDC_NO %in% GRDCtoremove_all), .N]

        ### Check changes in GSIM discharge data availability and flow regime over time ####
        GSIMflags <- rbindlist(list(GSIMtoremove_irartifacts,
                                    GSIMtoremove_winterIR,
                                    GSIMtoremove_unstableIR
        ))

        GSIMstatsdt_clean <- GSIMstatsdt[!(gsim_no %in%  GSIMflags[flag=='removed', gsim_no]),]


        mvars <- c('intermittent_o1800',
                   'intermittent_o1961',
                   'intermittent_o1971')
        alluv_formatGSIM <- melt(GSIMstatsdt_clean,
                                 id.vars = c('gsim_no',
                                             paste0('totalYears_kept_o',
                                                    c(1800, 1961, 1971))),
                                 measure.vars = mvars) %>%
            .[totalYears_kept_o1800 < 10 & variable %in% mvars, value := NA] %>%
            .[totalYears_kept_o1961 < 10 & variable %in% mvars[2:3], value := NA] %>%
            .[totalYears_kept_o1971 < 10 & variable %in% mvars[3], value := NA] %>%
            .[, count := .N, by=.(variable, value)]


        ### Check changes in GRDC discharge data availability and flow regime over time ####
        GRDCstatsdt_clean <- GRDCstatsdt[!(GRDC_NO %in% GRDCtoremove_all),]

        alluv_formatGRDC <- melt(GRDCstatsdt_clean,
                                 id.vars = c('GRDC_NO',
                                             paste0('totalYears_kept_o',
                                                    c(1800,1961, 1971))),
                                 measure.vars = mvars) %>%
            .[totalYears_kept_o1800 < 10 & variable %in% mvars, value := NA] %>%
            .[totalYears_kept_o1961 < 10 & variable %in% mvars[2:3], value := NA] %>%
            .[totalYears_kept_o1971 < 10 & variable %in% mvars[3], value := NA] %>%
            .[, count := .N, by=.(variable, value)]


        ###Analyze change in number of gauges with different intermittency criterion
        irsensi_format <- melt(rbind(GRDCstatsdt_clean, GSIMstatsdt_clean,
                                     use.names=TRUE, fill=T)[totalYears_kept_o1961 >= 10,],
                               id.vars = c('GRDC_NO', 'gsim_no'),
                               measure.vars = paste0('mDur_o', c(1800, 1961, 1971))) %>%
            .[!is.na(value) & value >0,] %>%
            setorder(variable, -value) %>%
            .[, cumcount := seq(.N), by=.(variable, is.na(GRDC_NO))]

        ggirsensi <- ggplot(irsensi_format, aes(x=value, y=cumcount,
                                                color=variable, linetype=is.na(GRDC_NO))) +
            geom_line(size=1.1) +
            coord_cartesian(expand=0, clip='off') +
            scale_x_sqrt(breaks=c(1, 5, 10, 30, 90, 180, 365),
                         labels=c(1, 5, 10, 30, 90, 180, 365)) +
            geom_vline(xintercept=c(1, 5)) +
            annotate(geom='text', x=c(1.7,6.5), y=150, angle=90,
                     label=c(sum(irsensi_format[value==1 & variable=='mDur_o1800',
                                                max(cumcount), by=is.na(GRDC_NO)]$V1),
                             sum(irsensi_format[value==5 & variable=='mDur_o1800',
                                                max(cumcount), by=is.na(GRDC_NO)]$V1))) +
            theme_classic()

        # plots <- grid.arrange(
        #     ggalluvium_gaugecount(dtformat = alluv_formatGRDC, alluvvar = 'GRDC_NO'),
        #     ggalluvium_gaugecount(dtformat = alluv_formatGSIM, alluvvar = 'gsim_no'),
        #     ggirsensi
        # )

        ### Bind GRDC and GSIM records ####################################
        databound <- rbind(GRDCstatsdt_clean,
                           GSIMstatsdt_clean,
                           use.names=TRUE, fill=T)

        return(list(#plots=plots,
                    data=databound,
                    flags=rbind(GRDCflags, GSIMflags,
                                use.names=TRUE, fill=T)))
    }

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
        .[dor_pc_pva>=5000, `:=`(flag = 'removed',
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
                                  "excluded from modeling  - perennial or undertermined",
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
