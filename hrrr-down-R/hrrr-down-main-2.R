#
#
#
#
# # # # # # # # #

####
## 1. INSPECT STATION VALUES FOR TIMESERIES 
# ATWATER 
ath_path = '/uufs/chpc.utah.edu/common/home/skiles-group4/otto/Atwater_MODIS_BBalbedo_comparison/Atwater_FullData_2022.csv'
df0 = snotel.data.ath.sub(ath_path, start_date = "2022-03-01", end_date = "2022-04-01")
df0 %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
# IRWIN
df0 = snotel.data.ath.sub(irwin_snotel21, start_date = "2021-03-01", end_date = "2021-05-01")
df0 %>% 
  mutate(station = Incoming_Solar_Wm2_1_Avg) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
# USGS IRWIN
ir_path_usgs = '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/irwin/Irwin_USGS.csv'
df0 = snotel.data.irwinusgs.sub(ir_path_usgs, start_date = "2022-03-26", end_date = "2022-04-01")
df0 %>% 
  mutate(station = visin) %>%
  ggplot(aes(x = loc_time)) + xlab("") +
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  scale_x_datetime(labels = scales::date_format("%b %d"), breaks = scales::date_breaks("1 day")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
# GRAND MESA
grand_path = '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/grand-mesa/GMSPC_wy2021.csv'
# check these values out!








####
## 2. REPRESENTATIVE DAYS
# # # # # # # # #
### FULL DAY 
# no save - 50 m - APRIL 1, 2021
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
single_date = "2021-04-01"
single_date = "2021-03-25"
single_date = "2021-03-30"
single_date = "2021-03-01"
single_date = "2021-02-18"
single_date = "2021-01-21"
hrrr.day.wrapper(day_string = single_date, hrrr_base, dpath, which.hrrr = 1, save_path = NULL, save_day = NULL)
jdf = snotel.merge.days(start_date = single_date, end_date = single_date, snotel_path = irwin_snotel21, multiday.vars=FALSE)
plot.day.station(jdf, span_p = 0.2)
plot.multiday.station(jdf)







# # # # # # # # #
### MULTIDAY RUN
# no save = SHORT TIME PERIOD ONLY
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
start_day = "2021-03-30"; end_day = "2021-04-01"
hrrr.multiday.wrapper(start_day = start_day, end_day = end_day, hrrr_base, dpath, which.hrrr = 1, 
                      save_path = NULL, save_day = NULL, save_time = NULL)
jdf = snotel.merge.days(start_date = start_day, end_date = end_day, snotel_path = irwin_snotel21, multiday.vars=TRUE)
plot.multiday.station(jdf)
plot.multiday.diff(jdf)






# # # # # # # # #
### MULTIDAY COMPILER FROM SAVED
# IRWIN
#
list.files(file.path(save_path, "hrrr-irwin-subset-day-21-bil"))
length(list.files(file.path(save_path, "hrrr-irwin-subset-day-21-bil")))/5
day_saved <- file.path(save_path, "hrrr-irwin-subset-day-21-bil") # IRWIN
dsf_irwin = save.merge.extract.five.sites(day_saved,  stat_pt = pts)
head(dsf_irwin)

# GRAND MESA
list.files(file.path(save_path, "hrrr-grandm-subset-day"))
day_saved <- file.path(save_path, "hrrr-grandm-subset-day") # IRWIN
dsf_gm = save.merge.extract.five.sites(day_saved,  stat_pt = pts_gm)
head(dsf_gm)

# stats
dsf_irwin %>% 
  mutate(month = format(as.POSIXct(time), "%m"),
         day = format(as.POSIXct(time), "%d"),) %>% 
  mutate(dswrf_diff = dswrf - Irwin,
         dswrf_z_diff = dswrf_z - Irwin,
         k_mod_diff = k_mod - Irwin,
         kn_mod_diff = kn_mod - Irwin,
         clear_diff = clear - Irwin) %>% 
  group_by(month) %>% 
  summarise(dswrf_mean = mean(dswrf_diff, na.rm=T),
            dswrfz_mean = mean(dswrf_z_diff, na.rm=T),
            k_mean = mean(k_mod_diff, na.rm=T),
            kn_mean = mean(kn_mod_diff, na.rm=T),
            clear_mean = mean(clear_diff, na.rm=T))


# stats
dsf_irwin %>% 
  mutate(month = format(as.POSIXct(time), "%m"),
         day = format(as.POSIXct(time), "%d"),) %>% 
  filter(Irwin<=0) %>% 
  mutate(dswrf_diff = dswrf - Irwin,
         dswrf_z_diff = dswrf_z - Irwin,
         k_mod_diff = k_mod - Irwin,
         kn_mod_diff = kn_mod - Irwin,
         clear_diff = clear - Irwin) %>% 
  group_by(month) %>% 
  summarise(dswrf_mean = mean(dswrf_diff, na.rm=T),
            dswrfz_mean = mean(dswrf_z_diff, na.rm=T),
            k_mean = mean(k_mod_diff, na.rm=T),
            kn_mean = mean(kn_mod_diff, na.rm=T),
            clear_mean = mean(clear_diff, na.rm=T))



### Overall RMSE - DSWRF
dsf_irwin %>% select(Irwin, dswrf) %>% na.omit() %>% 
  filter(Irwin>0.0) %>% #group_by(Irwin) %>% 
  summarise(rmse = sqrt(mean((Irwin - dswrf)^2)))
# Irwin RMSE: 406.2725
# kn_mod, clear
dsf_irwin %>% select(Irwin, kn_mod) %>% na.omit() %>% 
  filter(Irwin>0.0) %>% #group_by(Irwin) %>% 
  summarise(rmse = sqrt(mean((Irwin - kn_mod)^2)))
# 434.8
dsf_irwin %>% select(Irwin, clear) %>% na.omit() %>% 
  filter(Irwin>0.0) %>% #group_by(Irwin) %>% 
  summarise(rmse = sqrt(mean((Irwin - clear)^2)))
# 414.7



# PLOTS

# shapes_ = c("Station" = 17, 
#             "HRRR DSWRF" = 3,
#             "HRRR DSWRFz" = 1,
#             "EKDS model" = 1, 
#             "Clearsky" = "")
# plot time
colors = c("Station" = "forestgreen", 
           "HRRR DSWRF" = "grey20",
           "EKDS model" = "tomato3", 
           "Clearsky" = "deepskyblue1")
ltys = c("Station" = "solid", 
           "HRRR DSWRF" = "dotted",
           "EKDS model" = "solid", 
           "Clearsky" = "solid")
# PLOTS
df_table = dsf_irwin
# df_table = dsf_gm.     ## 3-18-20 // 2-10-14
start_date <- as.POSIXct("2021-02-12 00:00:00", tz='MST') # 3-18
end_date <- as.POSIXct("2021-02-19 23:59:59", tz = 'MST') # 3-20

df_table %>% 
  filter(time >= start_date & time <= end_date) %>% 
  mutate(station = Irwin) %>% ggplot(aes(x = time2)) + 
  # Clearsky background
  geom_line(aes(y = clear, colour = "Clearsky"), 
            size = 1.5, alpha = 0.2) +
  # Station line plot
  geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
  # HRRR DSWRF 
  geom_point(aes(y = dswrf, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.4, shape=3) +
  geom_line(aes(y = dswrf, colour = "HRRR DSWRF"), size = 0.75, alpha = 0.6, 
            linetype='dotted') +
  # EKDS model 
  geom_point(aes(y = kn_mod, colour = "EKDS model"), show.legend = TRUE, size = 2, alpha = 0.2) +
  geom_line(aes(y = kn_mod, colour = "EKDS model"), size = 0.5, alpha = 0.6) +
 
  # Color and formatting
  scale_color_manual(name = "", values = colors) +
  scale_x_datetime(labels = scales::date_format("%b %d"),
                   breaks = scales::date_breaks("1 day")) +
  xlab('') + ylab(expression(Irradiance ~ (W ~ m^-2))) + 
  theme_classic() + 
  theme( axis.text.x = element_text(angle = 50, hjust = 1, size=12),
         # legend.position = c(0.9, 0.9),
         plot.margin = margin(10, 20, 30, 30))


sv_name = "time-series-march1"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=10, height=5, units = "in", res=300)
p1
dev.off()
  


#

# dir.create("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/results/figs")














# different months! (good sample)