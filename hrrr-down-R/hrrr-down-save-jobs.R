#
#
# Save data - submit as a job!
#
#

# initialize
rm(list = ls());gc()
source('~/Documents/src/snowtools/hrrr_down_R/hrrr-down-functions-1.R', echo=FALSE)
strt = Sys.time()
# # # # # # # # #
### SAVE FILES !

# IRWIN - DONE file.path(save_path, "hrrr-irwin-subset-day"))
# stop()
#
print(dpath)
print(hrrr_base)
start_day = "2021-02-01"; end_day = "2021-06-30"
hrrr.multiday.save.wrapper(start_day = start_day, end_day = end_day,
                           hrrr_base, dpath, which.hrrr = 1,
                           save_path = NULL, save_day = save_path)
# stop()
# dir.create(file.path(save_path, "hrrr-irwin-subset-day-21-bil"))
# dir.exists(file.path(save_path, "hrrr-irwin-subset-day"))
# list.files(save_path)
# list.files(file.path(save_path, "hrrr-day-2"))
# list.files(file.path(save_path, "hrrr-irwin-subset-day-21-2"))
# list.files(file.path(save_path, "hrrr-atwater-subset-day"))

# ATWATER - RUNNING
# dir.create(file.path(save_path, "hrrr-atwater-subset-day-21"))
# dpath = '/uufs/chpc.utah.edu/common/home/skiles-group2/otto/GSLB_iSnobal_data/Outputs/Jordan/topo.nc'
# hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/'
# start_day = "2021-02-01"; end_day = "2021-05-01" # 2021 does not exist
# hrrr.multiday.save.wrapper(start_day = start_day, end_day = end_day,
#                            hrrr_base, dpath, which.hrrr = 2,
#                            save_path = NULL, save_day = save_path)
# #
# list.files(file.path(save_path, "hrrr-atwater-subset-day-21"))



# Grand Mesa
# # dir.create(file.path(save_path, "hrrr-grandm-subset-day"))
# dpath = dpath_gm
# start_day = "2021-02-01"; end_day = "2021-05-01"
# hrrr.multiday.save.wrapper(start_day = start_day, end_day = end_day,
#                            hrrr_base, dpath, which.hrrr = 1,
#                            save_path = NULL, save_day = save_path)
# 
# list.files(file.path(save_path, "hrrr-grandm-subset-day"))

# list.files(file.path(save_path, "hrrr-grandm-subset-day"))


# end time
Sys.time() - strt 