#
# Downscaling HRRR Downward Shortwave Radiation Elevational Downscaling
# R 4.03 Geospatial Packages (CHPC)
# HRRR v04
# MAIN - tests
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)

## TEST 1
# test individual
# April 1st 2022 Sunrise 6:52 AM MT (13:52 UTC) - Sunset 7:32 PM (19:32 MT or 2:32 UTC +1 day) 
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/hrrr.t13z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210328/hrrr.t07z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210328/hrrr.t10z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2021/hrrr.20210328/hrrr.t10z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/hrrr.20220401/hrrr.t14z.wrfsfcf06.grib2'
topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
read.hrrr(grib_path, dpoly, which.hrrr=1)
sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
clear.sky.sw(dem, dpoly, dtime, savepath=NULL) # returns sw3
hrrr.sw.topo(hsw, cld, dem, savepath=NULL) # hsw3
snotel.loc(irwin_path)
# original (joe) - new model
# negative values are reducing bias
plot( hsw03 - (dswrf * cos_inc) )
plot( hsw3 - (dswrf * cos_inc) )
extract(stack(hsw03 - (dswrf * cos_inc), hsw3 - (dswrf * cos_inc)), pts)
# SUCCESS!

## TEST 2
# try new scene
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/hrrr.t15z.wrfsfcf06.grib2'
hrrr.moment.wrapper(grib_path, dpath30, which.hrrr = 1, save_path = NULL)
# success

## TEST 3
# full day - no save (2 mins)
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
hrrr.day.wrapper(day_string = "2021-04-01", hrrr_base, dpath, which.hrrr = 1, save_path = NULL, save_day = NULL)
# success

# see  points at Irwin
snotel.loc(irwin_path)
hr_loc = (7:19)
plot(as.numeric(extract(dswrf_stk, pts))~hr_loc , type='l', col='black', ylim=c(0,900), ylab = "Incoming irradiance (Wm-2)", xlab = "Hour (local)")
lines(as.numeric(extract(hsw3_stk, pts))~hr_loc, col='blue')
lines(as.numeric(extract(hsw03_stk, pts))~hr_loc, col='green')
lines(as.numeric(extract(sw3_stk, pts))~hr_loc, col = 'firebrick', lty=2)

## TEST 4
# use higher resolution (3 mins)
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
hrrr.day.wrapper(day_string = "2021-04-01", hrrr_base, dpath30, which.hrrr = 1, save_path = NULL, save_day = NULL)
# SUCCESS

## TEST 5
# multiple days?


# compare to model
# DSWRF, illumination_angle, albedo, net_solar
hrrrjoe_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/erw_isnobal/wy2021/erw_hrrr_solar/run20210401/net_solar.nc'
net_solar = stack(hrrrjoe_path, varname="net_solar")
joe_dswrf = stack(hrrrjoe_path, varname="DSWRF")
alb = stack(hrrrjoe_path, varname="albedo")
crs(net_solar) <- "EPSG:32612";crs(alb) <- "EPSG:32612";crs(joe_dswrf) <- "EPSG:32612"

dsw0 <- net_solar / (1 - alb)
hr_full = 1:24
hr_loc = hr_full+6 # 6 hrs summer 7 winter
topo.bounds(dpath) 
snotel.loc(irwin_path)
plot(as.numeric(extract(dsw0, pts))~hr_full , type='l', col='black', ylim=c(0,1000), 
     ylab = "Incoming irradiance (Wm-2)", xlab = "Hour (local)")


lines(as.numeric(extract(hsw3_stk, pts))~hr_utc, col='blue')
lines(as.numeric(extract(dswrf_stk, pts))~hr_utc, col='red')
extract(dsw0, pts)[14:24] - extract(dswrf_stk, pts)[1:11]
extract(dsw0, pts)[14:24] - extract(hsw03_stk, pts)[1:11]


# SINGLE day snotel
# compare with snotel
# snotel_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2022/Irwin/2022-Irwin.csv'
snotel_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2021/Irwin/2021-Irwin.csv'
sdf = snotel.data.sub(snotel_path, start_date = "2021-04-01", end_date = "2021-04-01")
head(sdf)
colors = c("visin + nirin" = "forestgreen", 
           "visin*0.67 + nirin0.33" = "tomato3", 
           "visin" = "grey60",
           "nirin" = "grey60")
sdf %>% mutate(dswrf = visin + nirin) %>%
  mutate(dswrf2 = visin*0.67 + nirin*0.33) %>%
  ggplot(aes(x = loc_time)) + geom_line(aes( y = dswrf, color="visin + nirin"), linetype = 'solid') +
  geom_line(aes( y = dswrf2, color="visin*0.67 + nirin0.33"), linetype='solid') +
  geom_line(aes( y = visin, color="visin"), linetype='dashed') +
  geom_line(aes( y = nirin, color="nirin"),  linetype='dashed') +
  scale_x_datetime(labels = scales::date_format("%Y-%m-%d")) +
  xlab("") + ylab(expression(Irradiance ~ (W ~ m^-2))) +  theme_classic() +
  ggtitle("Irwin Plot") +
  scale_color_manual(name = "Net SW Down", values =  colors) +#c("deepskyblue", "forestgreen", "grey60", "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# create vars for 
ltime01 = gsub("\\D", "\\1", unlist(lapply(strsplit(names(dswrf_stk), "_"), function(x) x[[1]][1])))
start_date = "2021-04-01"
ltime0 = as.POSIXct(paste0(start_date, " ", ltime01, ":00:00") )
snotel.loc(irwin_path)
snow_date_df <- data.frame( loc_time = ltime0,
                            save_time = dtime_day,
                            zenith = zenith_day, 
            dswrf_inc = as.numeric(extract(dswrf_stk, pts)),
            dswrf_z = as.numeric(extract(dswrf_z_stk, pts)),
            kn_mod_inc = as.numeric(extract(hsw3_stk, pts)),
            kn_mod_z = as.numeric(extract(hsw3_z_stk, pts)),
            k_mod = as.numeric(extract(hsw03_stk, pts)),
            k_mod_z = as.numeric(extract(hsw03_z_stk, pts)),
            clearsky = as.numeric(extract(sw3_stk, pts)),
            clearsky_z = as.numeric(extract(sw3_stk_z, pts)))

snow_date_df[is.na(snow_date_df)] <- 0
head(snow_date_df)
# full join
joindf <- dplyr::left_join(sdf, snow_date_df, by='loc_time') 

# FINAL
start_date = "2021-04-01"
colors = c("Station_vis.67+nir.33" = "forestgreen", 
           "Station_vis" = "springgreen3", 
           "EKDS model" = "tomato3", 
           # "HRRR DSWRF" = "grey60",
           "Clearsky" = "deepskyblue1")
joindf %>% 
  mutate(station = visin * 0.67 + nirin * 0.33) %>%
  mutate(station2 = visin) %>%
  # mutate(station = visin + nirin) %>%
  ggplot(aes(x = loc_time)) + 
  # Line plot for 'Station' with correct color mapping
  geom_line(aes(y = station, colour = "Station_vis.67+nir.33"), size = 1.5, alpha = 0.4) +
  geom_line(aes(y = station2, colour = "Station_vis"), size = 1.5, alpha = 0.4) +
  
  # geom_point(aes(y = dswrf_z, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.2) +
  # geom_smooth(aes(y = dswrf_z, colour = "HRRR DSWRF"), se = FALSE, span = 0.5, size = 0.5, alpha = 0.6) +
  
  # Point and smoothing layers for 'kn_mod', 'k_mod', and 'clearsky' with specific colors
  geom_point(aes(y = kn_mod_z, colour = "EKDS model"), show.legend = TRUE, size = 2, alpha = 0.2) +
  geom_smooth(aes(y = kn_mod_z, colour = "EKDS model"), se = FALSE, span = 0.5, size = 0.5, alpha = 0.6) +
  
  geom_point(aes(y = clearsky, colour = "Clearsky"), show.legend = TRUE, size = 2, alpha = 0.2) +
  geom_smooth(aes(y = clearsky, colour = "Clearsky"), se = FALSE, span = 0.5, size = 0.5) +
  
  # # Manually setting colors in the legend
  # scale_color_manual(
  #   values = c("Station" = "forestgreen", 
  #              "k_mod" = "tomato3", 
  #              "Clearsky" = "deepskyblue1"),
  #   labels = c("Station", "EKS model", "Clearsky")
  # ) +
  scale_color_manual(name = "", values = colors) +
  
  # Formatting x-axis to display date and hour
  scale_x_datetime(labels = scales::date_format("%H:%M"), 
                   breaks = scales::date_breaks("5 hours")) +
  
  xlab(start_date) + 
  ylab(expression(Irradiance ~ (W ~ m^-2))) + 
  
  theme_classic() + 
  theme(legend.position = c(0.9, 0.9),  # Position legend inside the plot area
        plot.margin = margin(10, 20, 40, 40))  # Adjust plot margins


# joindf = read.csv('/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/results/join_df_20210401.csv')
# joindf$loc_time <- as.POSIXct(joindf$loc_time)
# write.csv(joindf, '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs/join_df_20210401.csv', 
          # row.names = F)


# # # # # # # # #
### FULL DAY 
# no save - 50 m - APRIL 1, 2021
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
single_date = "2021-04-01"
single_date = "2021-03-25"
single_date = "2021-03-30"
single_date = "2021-02-18"
single_date = "2021-01-21"
hrrr.day.wrapper(day_string = single_date, hrrr_base, dpath, which.hrrr = 1, save_path = NULL, save_day = NULL)
jdf = snotel.merge.days(start_date = single_date, end_date = single_date, snotel_path = irwin_snotel21, multiday.vars=FALSE)
plot.day.station(jdf, span_p = 0.5)
plot.day.station2(jdf)
plot.multiday.station(jdf)

# or original HRRR
hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2021' 
hrrr.day.wrapper(day_string = single_date, hrrr_base, dpath, which.hrrr = 2, save_path = NULL, save_day = NULL)
jdf = snotel.merge.days(start_date = single_date, end_date = single_date, snotel_path = irwin_snotel21, multiday.vars=FALSE)
plot.day.station(jdf)



# worked - issues with jdf table
# check variables
for (i in 1:12){plot(dswrf_stk[[i]], main = names(dswrf_stk[[i]]), col=heat.colors(11, 0.5), breaks = seq(0, 1000, 100));plot(pts, add=T,pch=1)}
for (i in 1:12){plot(hsw3_stk[[i]], main = names(hsw3_stk[[i]]), col=heat.colors(11, 0.5), breaks = seq(0, 1000, 100));plot(pts, add=T,pch=1)}
for (i in 1:12){plot(sw3_stk[[i]], main = names(sw3_stk[[i]]), col=heat.colors(11, 0.5), breaks = seq(0, 1000, 100));plot(pts, add=T,pch=1)}

# clouds are omitted! "2021-03-25" @ 11AM "local11_utc17_z47.81"
plot(dswrf_stk[[4]], col=heat.colors(255), main=names(dswrf_stk[[4]])) # why does this have terrain? - should not
# vs
plot(hsw03_stk[[4]], col=heat.colors(255), main=names(hsw03_stk[[4]]))

# # # # # # # # #
### MULTIDAY RUN
# no save
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
start_day = "2021-03-30"; end_day = "2021-04-01"
hrrr.multiday.wrapper(start_day = start_day, end_day = end_day, hrrr_base, dpath, which.hrrr = 1, 
                                  save_path = NULL, save_day = NULL, save_time = NULL)
jdf = snotel.merge.days(start_date = start_day, end_date = end_day, snotel_path = irwin_snotel21, multiday.vars=TRUE)
plot.multiday.station(jdf)
plot.multiday.diff(jdf)
#
# write.csv(jdf, '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/results/jdf_20210328_20210401_d50.csv',
#           row.names = F)
jdf = read.csv('/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/results/jdf_20210328_20210401_d50.csv')


# plot single variable
colors = c("Station" = "forestgreen", 
           # "EKDS model" = "tomato3", 
           "HRRR DSWRF" = "grey20")
           # "Clearsky" = "deepskyblue1")
jdf %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + 
  # Station line plot
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  # HRRR DSWRF with smoothing (ignoring NA's)
  geom_point(aes(y = dswrf_z, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.2) +
  geom_line(data = subset(jdf, !is.na(dswrf)), aes(y = dswrf, colour = "HRRR DSWRF"), size = 0.5, alpha = 0.6) +
  # # EKDS model with smoothing (ignoring NA's)
  # geom_point(aes(y = kn_mod_z, colour = "EKDS model"), show.legend = TRUE, size = 2, alpha = 0.2) +
  # geom_line(data = subset(jdf, !is.na(kn_mod_z)), aes(y = kn_mod_z, colour = "EKDS model"), size = 0.5, alpha = 0.6) +
  # scale_x_datetime(labels = scales::date_format("%b %d"),
  #                  breaks = scales::date_breaks("1 day")) +
  scale_x_datetime(labels = scales::date_format("%H:%M"), 
                   breaks = scales::date_breaks("1 hours")) +
  scale_color_manual(name = "", values = colors) +
  xlab('') + 
  ylab(expression(Irradiance ~ (W ~ m^-2))) + 
  theme_classic() + 
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         # legend.position = c(0.9, 0.9),
         plot.margin = margin(10, 20, 30, 30))



#### WATCH!
grib_date_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210328/'
grib_date_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/'
grib_date_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2021/hrrr.20210328'
grib_date_path = "/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/hrrr.20220401/"
# grib_date_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2021/hrrr.20210401'
# grib_date_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_CBR/hrrr.20210401/' # only includes Utah?
map.plotday.gif(grib_date_path, plot_h_level = 1, which.hrrr = 1) # 1:h - 2:hsw - 3:hsw3

## SAVE
# save data
stop()
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
start_day = "2021-02-01"; end_day = "2021-04-01"
hrrr.multiday.save.wrapper(start_day = start_day, end_day = end_day, 
                           hrrr_base, dpath, which.hrrr = 1, 
                           save_path = NULL, save_day = save_path)
stop()

### ATWATER
# TEST FOR ATWATER
ath_path = '/uufs/chpc.utah.edu/common/home/skiles-group4/otto/Atwater_MODIS_BBalbedo_comparison/Atwater_2023_data.csv'
ath_path = '/uufs/chpc.utah.edu/common/home/skiles-group4/otto/Atwater_MODIS_BBalbedo_comparison/Atwater_FullData_2022.csv'
# run for atwater
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
dpath = '/uufs/chpc.utah.edu/common/home/skiles-group2/otto/GSLB_iSnobal_data/Outputs/Jordan/topo.nc'
# ?? should minimize DEM extent?
hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/'
single_date = "2022-04-01"
gslb = FALSE; change  
hrrr.day.wrapper(day_string = single_date, hrrr_base, dpath, which.hrrr = 2, save_path = NULL, save_day = NULL)
ath_path = '/uufs/chpc.utah.edu/common/home/skiles-group4/otto/Atwater_MODIS_BBalbedo_comparison/Atwater_FullData_2022.csv'
jdf = snotel.merge.days(start_date = single_date, end_date = single_date, snotel_path = ath_path, 
                        multiday.vars=FALSE, is_ATH=TRUE)
plot.day.station(jdf)

# ATWATER MULTIDAY (inspect)
df0 = snotel.data.ath.sub(ath_path, start_date = "2022-03-01", end_date = "2022-04-01")
df0 %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))

### MULTIDAY RUN - ATWATER
# no save
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
dpath = '/uufs/chpc.utah.edu/common/home/skiles-group2/otto/GSLB_iSnobal_data/Outputs/Jordan/topo.nc'
hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/'
start_day = "2022-03-28"; end_day = "2022-04-01"
hrrr.multiday.wrapper(start_day = start_day, end_day = end_day, hrrr_base, dpath, which.hrrr = 2, 
                      save_path = NULL, save_day = NULL, save_time = NULL)
ath_path = '/uufs/chpc.utah.edu/common/home/skiles-group4/otto/Atwater_MODIS_BBalbedo_comparison/Atwater_FullData_2022.csv'
jdf = snotel.merge.days(start_date = start_day, end_date = end_day, snotel_path = ath_path, 
                        multiday.vars=FALSE, is_ATH=TRUE)

plot.multiday.station(jdf)
plot.multiday.diff(jdf)


jdf %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + 
  # geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  geom_line(aes(y=dswrf)) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() 

# write.csv(jdf, '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/ath/results/jdf_ATH_20220328_20220401_d100.csv',
#           row.names = F)
# jdf = read.csv('/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/ath/results/jdf_ATH_20220328_20220401_d100.csv')


create.dir('/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/ath/results')

list.files('/uufs/chpc.utah.edu/common/home/skiles-group4/otto')



df0 = snotel.data.ath.sub(irwin_snotel21, start_date = "2021-03-01", end_date = "2021-04-01")
df0 %>%
  # mutate(loc_time = as.POSIXct(TIMESTAMP)) %>% 
  mutate(visin = Incoming_Solar_Wm2_1_Avg) %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))

read.csv(irwin_snotel21) %>% 
  mutate(loc_time = as.POSIXct(TIMESTAMP, "%m/%d/%Y %H:%M", tz = "UTC")) %>% 
  filter(loc_time > as.POSIXct("2021-03-01") & loc_time < as.POSIXct("2021-04-01")) %>% 
  mutate(visin = Incoming_Solar_Wm2_1_Avg) %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))

# should have high bias -
# then use higher resolution DEM