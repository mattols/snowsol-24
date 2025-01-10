#
# Downscaling HRRR Downward Shortwave Radiation Elevational Downscaling
# R 4.03 Geospatial Packages (CHPC)
# HRRR v04
# FUNCTIONS
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

require(insol);require(raster);require(rgdal);require(rgeos)
require(ncdf4); require(sp);require(dplyr);require(ggplot2)


# # # # # PATHS # # # # # # # # # # # # # # # # # # # # # # # # # # # #
dpath = "/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/topo.nc" # 50 m
dpath30 = "/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw/erw_topo.grd" # 30 meter with larger extent
dpath_jord = '/uufs/chpc.utah.edu/common/home/skiles-group2/otto/GSLB_iSnobal_data/Outputs/Jordan/topo.nc' # 100 m
dpath_gm ="/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/grand-mesa/dem/gm_topo2.grd"
# irwin
hrrr_base = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR'
# hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2022'
save_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/results/agu-24/'

# GSL
# hrrr_base = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/'
# '/uufs/chpc.utah.edu/common/home/skiles-group3/GSLB_iSnobal_Outputs/MODIS/jordan_100m_isnobal_solar_albedo/wy2023/jordan_basin_100m/runYYYYMMDD/'

# snotel paths
irwin_snotel22 = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2022/Irwin/2022-Irwin.csv'
irwin_snotel21 = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2021/Irwin/2021-Irwin.csv'
irwin_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/site-locations/irwin_snotel.json'
# '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/site-locations/snotel_sites_3x3.json'

#######################################################################
# # FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################################

###############################################
# 0. DEM-BOUNDS
# set initial area and conditions (based on DEM)
topo.bounds <- function(dpath, subset=TRUE, gslb = FALSE){
  #
  # returns
  # dem, viewf, dpoly
  # read in DEM
  dem <<- raster(dpath, varname = 'dem')
  
  if(basename(dpath)=="topo.nc"){
    viewf <<- raster(dpath, varname = 'sky_view_factor')
    # VARS in .nc sky view & slope
    # viewf <<- raster(dpath, varname = 'sky_view_factor')
    # crs(viewf) <<- "EPSG:32613";crs(dem) <<- "EPSG:32613"
    # slp <<- raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/topo.nc", varname = 'slope')
    # crs(slp) <<- "EPSG:32613"
  }else{
    #
    viewf <<- brick(dpath)[[4]]
  }
  if( is.na( crs(viewf) )){
    if(gslb){crs(viewf) <<- "EPSG:32612";crs(dem) <<- "EPSG:32612"}else{crs(viewf) <<- "EPSG:32613";crs(dem) <<- "EPSG:32613"}
    # 
  }
  
  if(subset){
    # resample 30 meter dem for full area
    snotel.loc(irwin_path)
    if(gslb){ 
      dem <<- crop(dem, pts_a_buf9);viewf <<- crop(viewf, pts_a_buf9)
    }else{ # assume Irwin for now
      dem <<- crop(dem, pts_buf9);viewf <<- crop(viewf, pts_buf9)
      # dem <<- crop(dem, pts_gm_buf9);viewf <<- crop(viewf, pts_gm_buf9) # GRAND MESA
    }
    
  }
  
  # BOUNDARY - for computational spatial trasformations
  dpoly <<- as(extent(dem), 'SpatialPolygons')
  crs(dpoly) <<- crs(dem)
}


###############################################
# 1. HRRR-OBJ-TIME
# extract grib variables and time information
read.hrrr <- function(grib_path, dpoly, which.hrrr = 1){
  # 
  # HRRR grib file info
  # returns h, dtime (UTC), jdate
  hrrr <-  stack(grib_path)
  if(which.hrrr==1){
    names(hrrr) <- c('gmp', 'tempC', 'rh', 'u-wind', 'v-wind', 'prec06hr','prec01hr', 'cld0','cld', 'dswrf','vbdsf', 'vddsf')
  }else{
    names(hrrr) <- c('gmp', 'tempC', 'rh', 'u-wind', 'v-wind', 'prec06hr','prec01hr', 'cld0','cld', 'dswrf')
  }
  # clip and crop
  dpoly_proj <- spTransform(dpoly, proj4string(hrrr) )
  # reprojected 
  hrc <- crop(hrrr, dpoly_proj) #, snap='out')
  # hrc = crop(hrrr, buffer(dpoly_proj, 3.5e3)) # + 1 grid cell
  # h = projectRaster(hrc, dem) # change from 3000 m to 50 m - resampling near-nb by default
  h <<- projectRaster(hrc, crs = crs(dpoly) )
  # h <<- crop(h, dpoly)
  #
  # timestamp
  grib_hr <- gsub(".*t([0-9]+)z.*", '\\1', basename(grib_path) )
  day_date <-  as.Date(strsplit(basename(dirname(grib_path)), "\\.")[[1]][2], "%Y%m%d")
  moment_date_char <- paste0(as.character(day_date), " ", grib_hr, ":00:00")
  jdate <<- as.numeric(format(day_date, "%j"))
  dtime <<- as.POSIXct(moment_date_char, tz = "UTC") # UTC +1 for grid
}

###############################################
# 2. DEM-SW-DOWNSCALE
# downscale dswrf based on elevational relationship
sw.grid.elevation.resample2 <- function(h, dem, method='bilinear'){
  # 
  # elevational downscale resample of Dh
  # returns hsw
  #
  # create dissagregated raster for Dh
  
  # hd <- raster(dem); hd[] <- 1
  # # check factor
  # ff <-  res(h)[1]/res(dem)[1]
  # hdis  <-  disaggregate(h[['dswrf']], fact = ff) # 3000/50 = 60
  # hdd <-  resample(hdis, hd, method=method) # will not need to interpolate much but matches grid
  
  # model sw~elv using linear relationship (good approximation Arias-Ruiz et al., 2010)
  dat2 <- data.frame(x = values(h[['gmp']]), y =values(h[['dswrf']])) %>% na.omit()
  model <- lm(y ~ x, dat2)
  intercept <- coef(model)[1]  # Intercept (b0)
  slope <- coef(model)[2]     # Slope (b1)
  # print(paste('intercept:', intercept))
  # print(paste('slope:', slope))
  # clamp between range for cloudiness>
  slope <- pmax(pmin(slope, 2), 0)
  intercept <- pmax(pmin(intercept, 1200), 0)
  # predicted_y <- intercept + slope * new_x
  hsw <<- intercept + slope * dem
}

sw.grid.elevation.resample <- function(h, dem, method='bilinear'){
  # 
  # new option
  # NEW ELEVATIONAL DEPENDENCY!
  # 1. elevational relationship  
  # model sw~elv using linear relationship (good approximation Arias-Ruiz et al., 2010)
  dat2 <- data.frame(x = values(h[['gmp']]), y =values(h[['dswrf']])) %>% na.omit()
  model <- lm(y ~ x, dat2)
  intercept <- coef(model)[1]  # Intercept (b0)
  slope <- coef(model)[2]     # Slope (b1)
  slope <- pmax(pmin(slope, 1), 0)
  intercept <- pmax(pmin(intercept, 1000), 0)
  #
  # 2. check slope: if positive, run sample 
  # ELV-DOWN
  # if(slope > 0 & slope < 1 &  summary(model)$r.square > 0.2 ){
  if(FALSE){
    # CALCULATE ELEVATIONAL DEPENDENCY
    hd <- raster(dem); hd[] <- 1
    ff <-  res(h)[1]/res(dem)[1]
    hdis  <-  disaggregate(h[['dswrf']], fact = ff) # 3000/50 = 60
    hdd <-  resample(hdis, hd, method='ngb') # will not need to interpolate much but matches grid
    hsw <<- intercept + slope * dem
    # Correct for outlying points that could be cloud
    # residuals_ = residuals(model)
    # sigma_hat =  sd(residuals_) # check for other outliers?
    # # Identify outliers: residuals > 3 standard deviations away from 0
    # outlier_threshold <- 3 * sigma_hat
    # outliers <- which(abs(residuals_) > outlier_threshold)
    # # (optional) Cook's distance
    # cooks_dist <- cooks.distance(model)
    # influential_points <- which(cooks_dist > 4/length(dat2$x))
    # This threshold for Cook's distance (4/n) is a common rule of thumb for identifying points that may be disproportionately influencing the model
    #
    # identify pixels
    
    # UNABLE TO subset individual pixels due to interpolation
    
    # plot(y~x, data=dat2, main = dtime) 
    # abline(model)
    # # Plot points with outliers highlighted
    # points(dat2$x[outliers], dat2$y[outliers], col='red')
    # points(dat2$x[influential_points], dat2$y[influential_points], col='green')

    
    
  }else{
    ### ELSE USE (likely cloud)
    hsw <<- resample(h[['dswrf']], dem, method='bilinear')
  }
  
}



###############################################
# 3. CLEARSKY
# model insol with grid info
clear.sky.sw <- function(dem, dpoly, dtime, savepath=NULL){
  # model sw
  # returns:
  # zenith, Ib0, Id0, sh0, cosi 
  #   s - 'dem', 'viewf', 'Ib0', 'Id0', 'cos_inc', 'shade', 'Dh', 'Ib_terrain', 'Id_terrain', 'dswtotal'
  #
  # get coordinates
  if(gsub(".*proj=([a-z]+)\\s.*", '\\1', crs(dpoly) ) != "longlat"){
    pdpoly <- spTransform(dpoly, "EPSG:4326") 
  }
  lat_lon <- rowMeans(pdpoly@bbox)[2:1]
  # ddproj = projectRaster(dem, crs = crs("EPSG:4326"))
  # lat_lon <- c(round((ddproj@extent@ymax + ddproj@extent@ymin)/2,5), round((ddproj@extent@xmax + ddproj@extent@xmin)/2,5))
  # sun position and vector
  jd <-  JD(dtime)
  # timezone
  zchng = dtime; attr(zchng, "tzone") <- "US/Mountain"
  tmz = as.numeric(format(zchng, "%H")) - as.numeric(format(dtime, "%H"))
  tmz = -6
  # sunvector
  sv <-  sunvector(jd,lat_lon[1],lat_lon[2],tmz) # -6 if after first sunday in march 
  sp1 <- sunpos(sv)
  # ADD CHECK FOR SUNDOWN VALUES ?
  # daylight hours (zenith <= 90)
  # sp <-  ifelse(sp1[2]>=90, 90, sp1[2]) # sp1[which(sp1[,2]<=90),] 
  # sv <- sv[which(sp1[,2]<=90),]
  # sv <- c(1,-0,-0)
  # specify
  zenith <<-  sp1[2]
  azim <-  sp1[1]
  azimuth_eq <-  azim -180
  # run shade
  d_mat <- as.matrix(dem)
  dl <-  res(dem)[1]
  sh <- doshade(d_mat, sv, dl=dl) 
  sh0 <<- make.raster(sh, dem) 
  # zenith and incident angles
  s_a <- slope.aspect(dem)
  cos_inc <<- cos.slope(zenith, azimuth_eq, aspect = s_a[[2]], slope = s_a[[1]])
  # insolation at hour
  if(!exists('visibility')){insol.params()}
  height <- array(dem)
  Idirdif <-  insolation(zenith,jd,height,visibility,RH,tempK,O3,alphag)
  Ib <-  matrix(Idirdif[,1],nrow=nrow(dem),ncol=ncol(dem), byrow = T)
  Id <-  matrix(Idirdif[,2],nrow=nrow(dem),ncol=ncol(dem), byrow=T)
  # layers
  Ib0 <<- make.raster(Ib, dem)
  Id0 <<- make.raster(Id, dem) 
  #
  # final output - TOPOGRAPHIC SOLAR MODELS
  sw0 <<- Ib0 + Id0 * cos(radians(zenith))# Dh
  sw1 <<- Ib0 * cos_inc * sh0 # direct beam component
  sw2 <<- Id0 * cos(radians(zenith)) * viewf # diffuse component
  sw3 <<- sw1 + sw2 # direct and diffuse
  sw3_z <<- (Ib0 * cos(radians(zenith)) * sh0) + (Id0 * cos(radians(zenith)) * viewf)
  # stack
  s <<- stack(dem, viewf, Ib0, Id0, cos_inc, sh0, sw0, sw1, sw2, sw3, sw3_z)
  names(s) <<- c('dem', 'viewf', 'Ib0', 'Id0', 'cos_inc', 'shade', 'Dh', 'Ib_terrain', 'Id_terrain', 'dswtotal', 'dswtotal_z')
  if(!is.null(savepath)){
    # name
    dmoment_char <- paste0(gsub(":", "", gsub("\\s","_", gsub('-','',as.character(as.Date(dtime))) )), '_', format(dtime, "%H"), "utc_")
    zen_char <- gsub("\\.", "-",as.character(round(zenith, 2)) )
    sname = paste0("clearsky_", dmoment_char, "sz", zen_char, '.grd')
    writeRaster(s, file.path(savepath, "clearsky", sname))
  }
}

###############################################
# 4. HRRR-TOPO-SW
# new value for dswrf
hrrr.sw.topo <- function(hsw, cld, dem, savepath=NULL){
  #
  # final downscaled dswrf
  #
  # modify clouds
  # k_modify_clouds # .... # for another day ??? - should likely happen earlier in code
  
  # MODIFY HSW for zenith
  hsw_ <<- hsw 
  hsw_lim <- hsw
  hsw_lim[values(hsw_lim) > values(sw0)] <- sw0[values(hsw_lim) > values(sw0)] # clear sky is max
  hsw <<- hsw_lim / cos(radians(zenith))
   
  # clamp(hsw, max(values(sw0), na.rm=T))
  
  # calculate k and new_k
  sw.splitter(hsw, dem, sw0, Id0, zenith, jdate)
  
  # for comparison
  dswrf <<- resample(h[['dswrf']], dem, method='ngb') #
  dswrf_z <<- dswrf * cos(radians(zenith))
    
  # final output k
  Ib_h <- hsw * (1 - k)
  Id_h <- hsw * k
  hsw00 <<- Ib_h + Id_h # Dh
  hsw01 <<- Ib_h * cos_inc * sh0 # direct beam component
  hsw02 <<- Id_h * cos(radians(zenith)) * viewf # diffuse component
  hsw03 <<- hsw01 + hsw02 # direct and diffuse
  hsw03_z <<- ( Ib_h * cos(radians(zenith)) * sh0 ) + ( Id_h * cos(radians(zenith)) * viewf  ) 
  hrr_s <<- stack(dswrf, dswrf_z, hsw, k, Ib_h, Id_h, hsw00, hsw01, hsw02, hsw03, hsw03_z)
  names(hrr_s) <<- c('dswrf', 'dswrf_z', 'hsw','k','Ib_h', 'Id_h', 'Dh', 'Ib_terrain', 'Id_terrain', 'dswtotal', 'dswtot_zen')
  # final output new_k
  Ib_h <- hsw * (1 - new_k)
  Id_h <- hsw * new_k
  hsw0 <<- Ib_h + Id_h # Dh
  hsw1 <<- Ib_h * cos_inc * sh0 # direct beam component
  hsw2 <<- Id_h * cos(radians(zenith)) * viewf # diffuse component
  hsw3 <<- hsw1 + hsw2 # direct and diffuse
  hsw3_z <<- ( Ib_h * cos(radians(zenith)) * sh0 ) + ( Id_h * cos(radians(zenith)) * viewf  ) 
  hrr_s2 <<- stack(dswrf, dswrf_z, hsw, new_k, Ib_h, Id_h, hsw0, hsw1, hsw2, hsw3, hsw3_z)
  names(hrr_s2) <<- c('dswrf', 'dswrf_z', 'hsw','new_k', 'Ib_h', 'Id_h', 'Dh', 'Ib_terrain', 'Id_terrain', 'dswtotal', 'dswtot_zen')
  
  if(!is.null(savepath)){
    # name
    dmoment_char <- paste0(gsub(":", "", gsub("\\s","_", gsub('-','',as.character(as.Date(dtime))) )), '_', format(dtime, "%H"), "utc_")
    zen_char <- gsub("\\.", "-",as.character(round(zenith, 2)) )
    sname = paste0("hrrr-down_kG2_", dmoment_char, "sz",zen_char, '.grd')
    writeRaster(hrr_s, file.path(savepath, "hrrr-down", sname))
    # save second
    sname = paste0("hrrr-down_kn_", dmoment_char, "sz",zen_char, '.grd')
    writeRaster(hrr_s2, file.path(savepath, "hrrr-down", sname))
  }
  
}


###############################################
# 5. HRRR-MOD-DOWN-WRAPPER
# wrap generation of SW down for single grib file
hrrr.moment.wrapper <- function(grib_path, dpath, which.hrrr = 1, save_path = NULL){
  # run all for single grib file
  topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
  read.hrrr(grib_path, dpoly, which.hrrr) # returns: h, dtime (UTC), jdate 
  # intervene if values are 0
  # quick fix - not for saving
  # if(cellStats(h[['dswrf']], "max") <= 0){
  #   print("0% DSWRF - for grib scene")
  #   mainvars = c('dswrf', 'dswrf_z', 'sw3', 'sw3_z', 'hsw03', 'hsw3', 'hsw3_z', 'hsw3_z')
  #   secondvars = c('zenith', )
  #   for (x in mainvars) assign(x,NA, envir = .GlobalEnv)
  # }else{}
  # rest
  sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
  clear.sky.sw(dem, dpoly, dtime, savepath=save_path) # returns sw3
  hrrr.sw.topo(hsw, cld, dem, savepath=save_path) # hsw3
}

###############################################
# 6. HRRR-DAY TIMEWRAP
# execute over single day
hrrr.day.wrapper <- function(day_string = "2021-04-01", hrrr_base, dpath, which.hrrr = 1, save_path = NULL, save_day = NULL){
  #
  # run hrrr downscale for a specified day
  # full MODEL - need to run 24 hrs (for other vars)
  #
  # start
  # day_string = single_date
  strt = Sys.time()
  # read in dpoly
  if(!exists('dpoly')){topo.bounds(dpath)} # returns: dem, viewf, slp, dpoly
  dtime0 <- as.POSIXct(paste0(day_string, " 00:00:00"), tz = "UTC")
  # get coordinates
  if(gsub(".*proj=([a-z]+)\\s.*", '\\1', crs(dpoly) ) != "longlat"){
    pdpoly <- spTransform(dpoly, "EPSG:4326") 
  }
  lat_lon <- rowMeans(pdpoly@bbox)[2:1]
  # SUNRISE AND SET
  # April 1st 2022 Sunrise 6:52 AM MT (13:52 UTC) - Sunset 7:32 PM (19:32 MT or 2:32 UTC +1 day) 
  # MST is UTC-07:00 | MDT is UTC-06:00)
  # MST to MDT at 2 am MST to 3 am MDT on the second Sunday in March and returns at 2 am MDT to 1 am MST on the first Sunday in November
  # March 13 2am is switch
  # close but not as accurate
  # srs2 <- daylength(lat_lon[1], lat_lon[2], JD(dtime), -6) # jdate # JD(as.Date(dtime))
  # better (UTC) & MT
  sunr <- maptools::sunriset(cbind(lat_lon[2], lat_lon[1]), dtime0, direction = "sunrise", POSIXct.out=TRUE)
  suns <- maptools::sunriset(cbind(lat_lon[2], lat_lon[1]), dtime0, direction = "sunset", POSIXct.out=TRUE)
  sr_mt = sunr$time; attr(sr_mt, "tzone") <- "US/Mountain" # .POSIXct(sunr$time, tz="US/Mountain")
  ss_mt = suns$time; attr(ss_mt, "tzone") <- "US/Mountain"
  # find nearest hour
  # format(round(sunr$time, units= "hours"), "%H")
  # sun_start = lubridate::ceiling_date(sunr$time, unit = 'hours')
  # sun_end = lubridate::floor_date(suns$time, unit = 'hours')
  # local_start = lubridate::ceiling_date(sr_mt, unit = 'hours')
  # local_end = lubridate::floor_date(ss_mt, unit = 'hours')
  #
  sun_start = lubridate::floor_date(sunr$time, unit = 'hours')
  sun_end = lubridate::ceiling_date(suns$time, unit = 'hours')
  local_start = lubridate::floor_date(sr_mt, unit = 'hours')
  local_end = lubridate::ceiling_date(ss_mt, unit = 'hours')
  
  print(paste(" UTC: ", sunr$time, " - ", suns$time))
  print(paste("MT: ", sr_mt, " - ", ss_mt))
  print(paste("  HRRR hours", format(local_start, '%H'), "to", format(local_end, "%H")))
  
  grib_seq = c(format(sun_start, "%H"):23, 0:format(sun_end, "%H"))
  local_seq = format(local_start, "%H"):format(local_end, "%H")
  
  # check timechange
  tmz = as.numeric(format(local_start, "%H")) - as.numeric(format(sun_start, "%H"))
  
  # TIME is off - fails at 20 UTC - will not start at 7...
  # if(which.hrrr==1){grib_seq = local_seq}
  if(abs(tmz)< 7){
    grib_seq = local_seq # this works after March because 6th forecast hour compensates the UTC time
  } else{
    grib_seq = local_seq
    grib_seq = grib_seq + 1
  }
    # only for MT
  
  # subtract for grib forecast hour
  # grib_seq -6
  
  
  # run all scenes
  day_folder = paste0('hrrr.', format(dtime0, '%Y%m%d'))
  day_path = file.path(hrrr_base, day_folder)
  for (i in 1:length(grib_seq)){
    print(paste("...processing scene", i, "of", length(grib_seq)))
    
    # # specify grib tile
    # if (grib_seq[i] < grib_seq[1]){
    #   day_next <- paste0('hrrr.', format((dtime0 + 1*60*60*24), '%Y%m%d'))
    #   grib_path <- file.path(dirname(day_path), day_next, paste0('hrrr.t0', grib_seq[i], 'z.wrfsfcf06.grib2'))
    # }else{
    #   grib_path <- file.path(day_path, paste0('hrrr.t', grib_seq[i], 'z.wrfsfcf06.grib2'))
    # }
    # specify grib tile
    if (grib_seq[i] < grib_seq[1]){
      day_next <- paste0('hrrr.', format((dtime0 + 1*60*60*24), '%Y%m%d'))
      grib_path <- file.path(dirname(day_path), day_next, paste0('hrrr.t0', grib_seq[i], 'z.wrfsfcf06.grib2'))
    }else if(grib_seq[i] < 10){
      grib_path <- file.path(day_path, paste0('hrrr.t0', grib_seq[i], 'z.wrfsfcf06.grib2'))
    } else{
      grib_path <- file.path(day_path, paste0('hrrr.t', grib_seq[i], 'z.wrfsfcf06.grib2'))
    }
    # print(basename(grib_path))
    # print(grib_path)
    dtime1 <- as.POSIXct(paste0(dtime0, " ", grib_seq[i], ":00:00"),tz = 'US/Mountain'  )
    
    # run functions
    hrrr.moment.wrapper(grib_path, dpath, which.hrrr = which.hrrr, save_path = save_path)
    # keep main outputs in new stack
    if(grib_seq[i] == grib_seq[1]){
      
      dswrf_stk <<- dswrf #* cos_inc
      dswrf_z_stk <<- dswrf_z
      sw3_stk <<- sw3
      sw3_stk_z <<- sw0 #sw3_z
      hsw03_stk <<- hsw03
      hsw3_stk <<- hsw3
      hsw3_z_stk <<- hsw3_z
      hsw03_z_stk <<- hsw3_z
      # others
      zenith_day <<- zenith
      dtime_day <<- dtime1 
    }else{
      dswrf_stk <<- stack(dswrf_stk, dswrf ) #* cos_inc)
      dswrf_z_stk <<- stack(dswrf_z_stk, dswrf_z)
      sw3_stk <<- stack(sw3_stk, sw3)
      sw3_stk_z <<- stack(sw3_stk_z, sw3_z)
      hsw03_stk <<- stack(hsw03_stk, hsw03)
      hsw3_stk <<- stack(hsw3_stk, hsw3)
      hsw3_z_stk <<- stack(hsw3_z_stk, hsw3_z)
      hsw03_z_stk <<- stack(hsw03_z_stk, hsw3_z)
      # others
      zenith_day <<- c(zenith_day, zenith)
      dtime_day <<- c(dtime_day, dtime1)
    }
    # remove main vars
    rm(list=setdiff(ls(), c('strt', 'dtime0', 'grib_seq', 'day_folder', 'day_path', 'i',
                            'dpath', 'save_path', 'save_day', 'which.hrrr',
                            'dswrf_stk', 'sw3_stk', 'hsw03_stk', 'hsw3_stk',
                            'dswrf_z_stk', 'hsw3_z_stk', 'hsw03_z_stk',
                            'zenith_day', 'dtime_day') ) )
  }
  # set names # ONLY POST MARCH!
  # lltz <- ifelse(grib_seq[1]>grib_seq, grib_seq + (23-6), grib_seq-6)
  zchng = dtime; attr(zchng, "tzone") <- "US/Mountain"
  tmz = as.numeric(format(zchng, "%H")) - as.numeric(format(dtime, "%H"))
  # tmz = -6
  if (abs(tmz)< 7){
    lltz <- grib_seq+abs(tmz) # save UTC time
    names(dswrf_stk) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2)) # paste0('utc', grib_seq, '_local', lltz)
    names(dswrf_z_stk) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2))
    names(sw3_stk) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2))
    names(sw3_stk_z) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2))
    names(hsw03_stk) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2))
    names(hsw3_stk) <<- paste0('local', grib_seq, '_utc', lltz, '_z', round(zenith_day, 2))
  }else{ # BEFORE MARCH
    lltz <- grib_seq+abs(tmz) # save UTC time
    names(dswrf_stk) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2))
    names(dswrf_z_stk) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2)) # paste0('utc', grib_seq, '_local', lltz)
    names(sw3_stk) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2))
    names(sw3_stk_z) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2))
    names(hsw03_stk) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2))
    names(hsw3_stk) <<- paste0('local', grib_seq-1, '_utc', lltz, '_z', round(zenith_day, 2))
  }
  
  # save?
  if(!is.null(save_day)){
    # name
    dmoment_char <- gsub(":", "", gsub("\\s","_", gsub('-','',as.character(dtime0)) ))
    sname = paste0("hrrr-day_", dmoment_char, "_")
    # subfolder =   "hrrr-atwater-subset-day-21" # # "hrrr-day-2"
    # subfolder = "hrrr-grandm-subset-day"
    subfolder = "hrrr-irwin-subset-day-21-bil"
    writeRaster(dswrf_stk, file.path(save_day, subfolder, paste0(sname, 'dswrf.grd') ) )
    writeRaster(dswrf_z_stk, file.path(save_day, subfolder, paste0(sname, 'dswrf_z.grd') ) )
    # writeRaster(sw3_stk, file.path(save_day, subfolder, paste0(sname, 'clearsky.grd') ) )
    writeRaster(sw3_stk_z, file.path(save_day, subfolder, paste0(sname, 'clearsky.grd') ) )
    writeRaster(hsw03_z_stk, file.path(save_day, subfolder, paste0(sname, 'hrrr_down_kG2.grd') ) )
    writeRaster(hsw3_z_stk, file.path(save_day, subfolder, paste0(sname, 'hrrr_down_kn.grd') ) )
  }
  # end time
  Sys.time() - strt 
  
}


###############################################
# 7. HRRR-MULTIDAY TIMEWRAP
# execute over single day - SHORT timeperiods only!
hrrr.multiday.wrapper <- function(start_day = "2021-03-27", end_day = "2021-04-01", hrrr_base, dpath, which.hrrr = 1, 
                                  save_path = NULL, save_day = NULL, save_time = NULL){
  #
  # specify start and end days as strings
  #
  # start
  strt = Sys.time()
  
  # dates list
  dates_ls <- seq(as.Date(start_day), as.Date(end_day), by="day")
  # check length for computation
  if((length(dates_ls)>5) ){stop("! Shorten time period to 8 days or less!")} # & !is.null(save_day)
  
  for(dl in 1:length(dates_ls)){
    # run current day
    current_day = as.character(dates_ls[dl])
    hrrr.day.wrapper(day_string = current_day, hrrr_base, dpath, which.hrrr = which.hrrr, save_path = save_path, save_day = save_day)
    
    # keep main outputs in new stack
    if(dates_ls[dl] == dates_ls[1]){
      # save stacks
      dswrf_stkm <<- dswrf_stk
      dswrf_z_stkm <<- dswrf_z_stk
      sw3_stkm <<- sw3_stk
      sw3_stk_zm <<- sw3_stk_z
      hsw03_stkm <<- hsw03_stk
      hsw3_stkm <<- hsw3_stk
      hsw3_z_stkm <<- hsw3_z_stk
      hsw03_z_stkm <<- hsw03_z_stk
      # others
      zenith_daym <<- zenith_day
      dtime_daym <<- dtime_day 
    }else{
      dswrf_stkm <<- stack(dswrf_stkm,dswrf_stk)
      dswrf_z_stkm <<- stack(dswrf_z_stkm, dswrf_z_stk)
      sw3_stkm <<- stack(sw3_stkm, sw3_stk)
      sw3_stk_zm <<- stack(sw3_stk_zm, sw3_stk_z)
      hsw03_stkm <<- stack(hsw03_stkm, hsw03_stk)
      hsw3_stkm <<- stack(hsw3_stkm, hsw3_stk)
      hsw3_z_stkm <<- stack(hsw3_z_stkm, hsw3_z_stk)
      hsw03_z_stkm <<- stack(hsw03_z_stkm, hsw03_z_stk)
      # others
      zenith_daym <<- c(zenith_daym, zenith_day)
      dtime_daym <<- c(dtime_daym, dtime_day)
    }
    # remove main vars
    rm(list=setdiff(ls(), c('dates_ls', 'dl', 'which.hrrr', 'strt',
                            'hrrr_base', 'dpath', 'save_path', 'save_day', 'save_time',
                            'dswrf_stkm', 'sw3_stkm', 'hsw03_stkm', 'hsw3_stkm',
                            'dswrf_z_stkm', 'hsw3_z_stkm', 'hsw03_z_stkm',
                            'zenith_daym', 'dtime_daym') ) )
  }
  
  print( paste("Completed downscale from:", start_day, "to", end_day) )
  # end time
  Sys.time() - strt 
}


###############################################
# 9. HRRR-DAY TIMEWRAP
# execute over several days
hrrr.multiday.save.wrapper <- function(start_day = "2021-03-01", end_day = "2021-04-01", hrrr_base, dpath, which.hrrr = 1, 
                                  save_path = NULL, save_day){
  #
  # specify start and end days as strings
  #
  # start
  strt = Sys.time()
  
  # dates list
  dates_ls <- seq(as.Date(start_day), as.Date(end_day), by="day")

  for(dl in 1:length(dates_ls)){
    # run current day
    current_day = as.character(dates_ls[dl])
    hrrr.day.wrapper(day_string = current_day, hrrr_base, dpath, 
                     which.hrrr = which.hrrr, save_path = save_path, save_day = save_day)
    
    # remove main vars
    rm(list=setdiff(ls(), c('dates_ls', 'dl', 'which.hrrr', 'strt',
                            'hrrr_base', 'dpath', 'save_path', 'save_day', 'save_time') ) )
  }
  
  print( paste("Completed downscale from:", start_day, "to", end_day) )
  # end time
  Sys.time() - strt 
}




###############################################
# 10. SNOTEL FUNCTIONS
# extract snotel location and timeseries data
snotel.loc <- function(snotel_path){
  #
  # should be json file
  # convert to sp spatial vector object
  #
  
  # create for atwater
  pts_a <<- SpatialPoints(coords = cbind(-111.63778, 40.59148), proj4string = CRS('EPSG:4326') )
  pts_a <<- spTransform(pts_a, proj4string(dem))
  pts_a_buf9 <<- raster::buffer(pts_a, 3.5e3) # 3500 m surrounding cells
  
  # create for Grand Mesa
  pts_gm <<- SpatialPoints(coords = cbind(-108.061435, 39.050802), 
                           proj4string = CRS('EPSG:4326') )
  pts_gm <<- spTransform(pts_gm, proj4string(dem))
  pts_gm_buf9 <<- raster::buffer(pts_gm, 3.5e3) # 3500 m surrounding cells
  
  # Sentaor Beck
  # SASP 
  sasp = cbind(37.906914, 107.711322) # https://snowstudies.org/swamp-angel-study-plot/
  # SBB 
  # cbind(37.907028, -107.711175) # less accurate
  sbb = cbind(37.906914, -107.726265) # https://snowstudies.org/senator-beck-study-plot/
  
  # extract lat lon
  con = file(snotel_path)
  line = readLines(con)
  lon0 = trimws(grep("lon", line, value = T))
  lat0 = trimws(grep("lat", line, value = T))
  x = as.numeric(strsplit(strsplit(lon0, ": ")[[1]][2], ',')[[1]][1])
  y = as.numeric(strsplit(strsplit(lat0, ": ")[[1]][2], ',')[[1]][1])
  # create points
  # pts <- cbind(x, y)
  pts <<- SpatialPoints(coords = cbind(x, y), proj4string = CRS(proj4string(dem)) )
  pts_buf9 <<- raster::buffer(pts, 3.5e3) # 3500 m surrounding cells
}

# snotel data.frame
snotel.data.sub <- function(snotel_path, start_date = "2022-03-26", end_date = "2022-04-01"){
  #
  # snotel_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2022/Irwin/2022-Irwin.csv'
  idf = read.csv(snotel_path, row.names = NULL)
  names(idf) <- c('time', 'temp', 'visin', 'nirin', 'visout', 'nirout')
  # max 1300 min 0 - no out when no in 
  idf$visin[idf$visin > 1300] = 1300; idf$visin[idf$visin < 0] = 0
  idf$nirin[idf$nirin > 1300] = 1300; idf$nirin[idf$nirin < 0] = 0
  idf$visout[idf$visin <= 0] = 0; idf$visin[idf$visin < 0] = 0
  idf$nirout[idf$nirin <= 0] = 0; idf$nirin[idf$nirin < 0] = 0
  # ## no OUTCOMING larger than INCOMING
  # irwin_pyra['VIS_OUT'] = irwin_pyra.apply(lambda r: r['VIS_OUT'] if r['VIS_OUT'] < r['VIS_IN'] else r['VIS_IN'], axis=1)
  # irwin_pyra['NIR_OUT'] = irwin_pyra.apply(lambda r: r['NIR_OUT'] if r['NIR_OUT'] < r['NIR_IN'] else r['NIR_IN'], axis=1)

  # as.POSIXct(idf$time[1], "%m/%d/%Y %H:%M", tz = "UTC")
  idf$time = as.POSIXct(idf$time, "%m/%d/%Y %H:%M", tz = "UTC")
  
  # Step 1: Define the start and end dates for the range
  # start_date <- as.POSIXct("2022-03-26 00:00:00") # 19 for two weeks
  # end_date <- as.POSIXct("2022-04-01 23:59:59")
  start_date <- as.POSIXct(paste0(start_date, "00:00:00") )
  end_date <- as.POSIXct(paste0(end_date, "23:59:59") )
  strt_utc = start_date; attr(strt_utc, "tzone") <- "UTC" # "US/Mountain"
  end_utc = end_date; attr(end_utc, "tzone") <- "UTC"
  
  # Step 2: Subset the data based on the date range
  idf2 <- idf[idf$time >= strt_utc & idf$time <= end_utc, ]
  
  idf2$loc_time <- idf2$time; attr(idf2$loc_time, "tzone") <- "US/Mountain"
  
  return(idf2)
  # idf2 %>% mutate(dswrf = visin + nirin) %>% 
  #   ggplot(aes(x = time, y = dswrf)) + geom_line() +
  #   theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  #   scale_x_datetime(labels = scales::date_format("%Y-%m-%d")) +
  #   xlab("")
}

# snotel data.frame
snotel.data.irwinusgs.sub <- function(snotel_path, start_date = "2022-03-26", end_date = "2022-04-01"){
  #
  # Data provided for site 385315107063001
  #            TS   parameter     Description
  #        307700       72175     Longwave radiation, downward intensity, watts per square meter
  #        307704       72186     Shortwave radiation, downward intensity, watts per square meter
  #        307705       72253     Soil temperature, degrees Celsius, [20 cm depth CS655]
  #        307706       00036     Wind direction, degrees clockwise from true north
  #        307708       74207     Moisture content, soil, volumetric, percent of total volume, [20 cm depth CS655]
  #        307711       72189     Snow depth, meters
  #        307712       72185     Shortwave radiation, upward intensity, watts per square meter
  #        307714       72393     Liquid water content, snowpack, percent of total volume, [10 cm above soil surface]
  #        307715       72341     Water content of snow, millimeters
  #        307716       72253     Soil temperature, degrees Celsius, [5 cm depth CS655]
  #        307717       00025     Barometric pressure, millimeters of mercury
  #        307719       72392     Snowpack density, kilograms per cubic meter, [10 cm above soil surface]
  #        307720       72174     Longwave radiation, upward intensity, watts per square meter
  #        307721       00052     Relative humidity, percent
  #        307722       72394     Mass flux density of drifting snow particles, grams per square meter per second
  #        307723       00020     Temperature, air, degrees Celsius
  #        307726       00035     Wind speed, miles per hour
  #        307727       74207     Moisture content, soil, volumetric, percent of total volume, [5 cm depth CS655]
  #        309043       72405     Surface temperature, non-contact, degrees Celsius

  idf = read.csv(snotel_path, row.names = NULL)
  idf = idf %>% select(datetime, X307700_72175, X307704_72186) %>% 
    rename(Longwave_down = X307700_72175  )   %>% 
    rename(Shortwave_down = X307704_72186) %>% 
    mutate(visin = Shortwave_down) %>% 
    mutate(lwin = Shortwave_down) %>% select(datetime,visin, lwin)
  # Soil_t20 = X307705_72253
  # Wind_dir = X307706_00036
  # Moisture_soil20 = X307708_74207
  # Snow_depth = X307711_72189
  # Shortwave_up = X307712_72185
  # Liquid_water_snow10 = X307714_72393
  # Water_content_snow = X307715_72341
  # Soil_t5 = X307716_72253
  # Barometric_pressure = X307717_00025
  # Snow_dense_10 = X307719_72392
  # Longwave_up = X307720_72174
  # RH = X307721_00052
  # Mass_flux_snow = x307722_72394
  # Tair = x307723_00020
  # Wind_speed = X307726_00035
  # Moisture_soil20 = x307727_74207
  # Surface_temp = x309043_72405

  idf$visin[idf$visin > 1300] = 1300; idf$visin[idf$visin < 0] = 0
  
  # as.POSIXct(idf$time[1], "%m/%d/%Y %H:%M", tz = "UTC")
  idf$time = as.POSIXct(idf$datetime, "%m/%d/%y %H:%M", tz = "MST")
  
  # Step 1: Define the start and end dates for the range
  # start_date <- as.POSIXct("2022-03-26 00:00:00") # 19 for two weeks
  # end_date <- as.POSIXct("2022-04-01 23:59:59")
  start_date <- as.POSIXct(paste0(start_date, "00:00:00"), tz = "MST" )
  end_date <- as.POSIXct(paste0(end_date, "23:59:59"), tz = "MST" )
  # strt_utc = start_date; attr(strt_utc, "tzone") <- "UTC" # "US/Mountain"
  # end_utc = end_date; attr(end_utc, "tzone") <- "UTC"
  
  # Step 2: Subset the data based on the date range
  idf2 <- idf[idf$time >= start_date & idf$time <= end_date, ]
  
  idf2$loc_time <- idf2$time #; attr(idf2$loc_time, "tzone") <- "US/Mountain"
  
  return(idf2)
  # idf2 %>% mutate(dswrf = visin + nirin) %>% 
  #   ggplot(aes(x = time, y = dswrf)) + geom_line() +
  #   theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  #   scale_x_datetime(labels = scales::date_format("%Y-%m-%d")) +
  #   xlab("")
}


# snotel data.frame
snotel.data.ath.sub <- function(snotel_path, start_date = "2022-03-26", end_date = "2022-04-01"){
  #
  # snotel_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2022/Irwin/2022-Irwin.csv'
  idf = read.csv(snotel_path)
  # idf$TIMESTAMP = as.POSIXct(idf$TIMESTAMP, "%m/%d/%Y %H:%M", tz = "UTC")
  idf$TIMESTAMP = as.POSIXct(idf$TIMESTAMP, "%m/%d/%Y %H:%M", tz = "US/Mountain")
  # start date
  start_date <- as.POSIXct(paste0(start_date, "00:00:00") )
  end_date <- as.POSIXct(paste0(end_date, "23:59:59") )
  # strt_utc = start_date; attr(strt_utc, "tzone") <- "UTC" # "US/Mountain"
  # end_utc = end_date; attr(end_utc, "tzone") <- "UTC"
  # idf2 <- idf[idf$TIMESTAMP >= strt_utc & idf$TIMESTAMP <= end_utc, ]
  idf2 <- idf[idf$TIMESTAMP >= start_date & idf$TIMESTAMP <= end_date, ]
  idf2$loc_time <- idf2$TIMESTAMP; #attr(idf2$loc_time, "tzone") <- "US/Mountain"
  idf2$visin = idf2$inVis
  return(idf2)
}

###############################################
# 11. HRRR-NETSOLAR - ORIGINAL FILES
# net_solar var for spatial comparison
snotel.merge.days <- function(start_date = "2021-04-01", end_date = "2021-04-01", snotel_path, 
                              multiday.vars=FALSE, is_ATH=FALSE){
  #
  # check 
  if(!any(c(exists('hsw03_stk') , exists('hsw03_stkm')))){stop("! Run model first !")}
  
  snotel.loc(irwin_path) # points
  
  if(multiday.vars){
    # save stacks
    dswrf_stk <- dswrf_stkm
    dswrf_z_stk <- dswrf_z_stkm
    sw3_stk <- sw3_stkm
    sw3_stk_z <- sw3_stk_zm
    hsw03_stk <- hsw03_stkm
    hsw3_stk <- hsw3_stkm
    hsw3_z_stk <- hsw3_z_stkm
    hsw03_z_stk <- hsw03_z_stkm
    # others
    zenith_day <- zenith_daym
    dtime_day <- dtime_daym
    # time
    if(is_ATH){
      # ndays = as.numeric((as.Date(end_date) - as.Date(start_date) ) +1) # number of days
      all_dates = seq(as.Date(start_date), as.Date(end_date), by = 'day')
      ndays = length(all_dates)
      ltime01 = gsub("\\D", "\\1", unlist(lapply(strsplit(names(dswrf_stkm), "_"), function(x) x[[1]][1])))
      ltime0 = as.POSIXct(paste0(rep(all_dates, ndays), " ", ltime01, ":00:00"))
    }else{ltime0 = dtime_daym}
    
  }else{
    ltime01 = gsub("\\D", "\\1", unlist(lapply(strsplit(names(dswrf_stk), "_"), function(x) x[[1]][1])))
    ltime0 = as.POSIXct(paste0(start_date, " ", ltime01, ":00:00") )
  }
  # snotel data
  # snotel_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2021/Irwin/2021-Irwin.csv'
  if(is_ATH){
    sdf = snotel.data.ath.sub(snotel_path, start_date = start_date, end_date = end_date)
    stat_pt = pts_a
  }else{
    sdf = snotel.data.sub(snotel_path, start_date = start_date, end_date = end_date)
    stat_pt = pts
  }
  
  
  # point data and extract by point
  
  snow_date_df <- data.frame( loc_time = ltime0,
                              save_time = dtime_day,
                              zenith = zenith_day, 
                              dswrf = as.numeric(extract(dswrf_stk, stat_pt)),      # Original HRRR with neareast neighbot
                              dswrf_z = as.numeric(extract(dswrf_z_stk, stat_pt)),  # HRRR on horizontal plane
                              kn_mod_inc = as.numeric(extract(hsw3_stk, stat_pt)),  # K-model (new) on inclined plane
                              kn_mod_z = as.numeric(extract(hsw3_z_stk, stat_pt)),  # K-model (new) on flat plane
                              k_mod = as.numeric(extract(hsw03_stk, stat_pt)),      # K-model (G2 local) on inclined
                              k_mod_z = as.numeric(extract(hsw03_z_stk, stat_pt)),  # K-model (G2 local) on flat plane
                              clearsky = as.numeric(extract(sw3_stk, stat_pt)),     # clearsky model
                              clearsky_z = as.numeric(extract(sw3_stk_z, stat_pt))) # clearsky on flat 

  snow_date_df[is.na(snow_date_df)] <- 0
  # head(snow_date_df)
  # full join
  joindf <- dplyr::left_join(sdf, snow_date_df, by='loc_time') 
  
  return(joindf)
}


## DF FROM SAVED
save.merge.extract.five.sites <- function(day_saved,  stat_pt){
  #
  # defaults for irwin points
  #
  # variable file lists
  dswrf_save <- list.files(day_saved, pattern = "dswrf.grd", full.names = T)
  dswrf_z_save <- list.files(day_saved, pattern = "dswrf_z.grd", full.names = T)
  k_save <- list.files(day_saved, pattern = "hrrr_down_kG2.grd", full.names = T)
  kn_save <- list.files(day_saved, pattern = "hrrr_down_kn.grd", full.names = T)
  clear_save <- list.files(day_saved, pattern = "clearsky.grd", full.names = T)
  
  # get DAY dates (from dswrf)
  # strsplit(list.files(day_saved, pattern = "dswrf.grd")[1], "_")
  save_dates <- unlist(lapply(dswrf_save, function(x) strsplit(x, "_")[[1]][2] ))
  
  # could create blank but maybe subset is better?
  # df00 = data.frame()
  dfive = read.csv("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/snotel-sites/snotel-five-sites-022021-062021.csv")
  split_time = unlist(lapply(1:length(dfive$time), function(x) strsplit(dfive$time[x], " M")[[1]][1]))
  split_tz = unlist(lapply(1:length(dfive$time), function(x) strsplit(dfive$time[x], " M")[[1]][2]))
  dfive$time = as.POSIXct(split_time, "%m/%d/%Y %H:%M", tz = "MST")
  # Adjust for MDT
  dfive$time2 = dfive$time
  dfive$time2[which(split_tz=="DT")] = dfive$time2[which(split_tz=="DT")] + (60*60)
  
  dfive$dirty_time = ifelse(split_tz=="ST", -1, 3) # -1 , 2
  # dfive$dirty_time = 0
  
  dfive$dswrf <- 0
  dfive$dswrf_z <- 0
  dfive$k_mod <- 0
  dfive$kn_mod <- 0
  dfive$clear <- 0
  dfive$zenith <- 0
  
  snotel.loc(irwin_path)
  
  # add more for other sites - irwin is now
  for (i in 1:length(dswrf_save)){
    # extract daily variables
    # pts    # Irwin
    # pts_a  # Atwater
    # pts_gm # Grand Mesa # CANNOT RUN
    # # test
    # i = 1
    # save_date = "2021-02-01"
    # i = 60
    # save_date = "2021-04-01"
    print(paste(i, "of", length(dswrf_save), " - ", save_dates[i]))
    
    # get_variables for current day
    dswrf_stk00    <- stack(dswrf_save[i]) 
    dswrf_z_stk00  <- stack(dswrf_z_save[i]) 
    k_stkmod00    <- stack(k_save[i]) 
    kn_mod00      <- stack(kn_save[i]) 
    clear00      <- stack(clear_save[i]) 
    
    # get datetime of images
    save_date = as.Date(save_dates[i], "%Y%m%d")
    ltime01 = gsub("\\D", "\\1", unlist(lapply(strsplit(names(dswrf_stk00), "_"), function(x) x[[1]][1])))
    # SHOULD MATCH? - needs adjustment???
    ltime0 = as.POSIXct(paste0(save_date, " ", ltime01, ":00:00"), tz = "US/Mountain" )
    
    zen0 = as.numeric(gsub("z", "", unlist(lapply(strsplit(names(dswrf_stk00), "_"), function(x) x[3]))))
    
    # range of values within final table
    # dfive[ match(ltime0, dfive$time), ]
    # IRWIN
    # stat_pt = pts
    range_pts = match(ltime0, dfive$time2) + dfive$dirty_time[match(ltime0, dfive$time2)]
    # range_pts = range_pts[!is.na(range_pts)]
    dfive$dswrf[ range_pts[!is.na(range_pts)] ]   <- as.numeric(extract(dswrf_stk00, stat_pt))[!is.na(range_pts)]
    dfive$dswrf_z[ range_pts[!is.na(range_pts)] ] <- as.numeric(extract(dswrf_z_stk00, stat_pt))[!is.na(range_pts)]
    dfive$k_mod[ range_pts[!is.na(range_pts)] ]   <- as.numeric(extract(k_stkmod00, stat_pt))[!is.na(range_pts)]
    dfive$kn_mod[range_pts[!is.na(range_pts)] ]  <- as.numeric(extract(kn_mod00, stat_pt))[!is.na(range_pts)]
    dfive$clear[range_pts[!is.na(range_pts)] ]  <- as.numeric(extract(clear00, stat_pt))[!is.na(range_pts)]
    dfive$zenith[range_pts[!is.na(range_pts)] ]  <- zen0[!is.na(range_pts)]
    
    # test
    # head(dfive, 20)[,6:9] # need to subtract an hour? MST
    # dfive[1415:1435,6:9] # need to add 2 hours? MDT
    # USES DIRT TIME for match???
    
    
  }
  
  return(dfive)
}


###############################################
# 12. 
#







###############################################
# PLOTTING FUNCTIONS # # # # # # # # # # # # #
###############################################

# plot single day
plot.day.station <- function(jdf, span_p = 0.5){
  
  # if(class(jdf$loc_time)[1]!='POSIXct'){
  #   jdf$loc_time <- as.POSIXct(jdf$loc_time)
  # }
  
  colors = c("Station" = "forestgreen", 
             "EKDS model" = "tomato3", 
             "HRRR DSWRF" = "grey60",
             "Clearsky" = "deepskyblue1")
  jdf %>% 
    # mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
    # mutate(station = visin * 0.67 + nirin * 0.33) %>%
    mutate(station = visin) %>%
    # mutate(station = visin + nirin) %>%
    ggplot(aes(x = loc_time)) + 
    # Line plot for 'Station' with correct color mapping
    geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
    
    geom_point(aes(y = dswrf, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_smooth(aes(y = dswrf, colour = "HRRR DSWRF"), se = FALSE, span = span_p, size = 0.5, alpha = 0.6) +
    
    # Point and smoothing layers for 'kn_mod', 'k_mod', and 'clearsky' with specific colors
    geom_point(aes(y = kn_mod_z, colour = "EKDS model"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_smooth(aes(y = kn_mod_z, colour = "EKDS model"), se = FALSE, span = span_p, size = 0.5, alpha = 0.6) +
    
    geom_point(aes(y = clearsky, colour = "Clearsky"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_smooth(aes(y = clearsky, colour = "Clearsky"), se = FALSE, span = span_p, size = 0.5) +
    scale_color_manual(name = "", values = colors) +
    scale_x_datetime(labels = scales::date_format("%H:%M"), 
                     breaks = scales::date_breaks("5 hours"), timezone = "US/Mountain") +
    xlab('') + 
    ylab(expression(Irradiance ~ (W ~ m^-2))) + 
    theme_classic() + 
    theme(legend.position = c(0.9, 0.9),  # Position legend inside the plot area
          plot.margin = margin(10, 20, 30, 30))  # Adjust plot margins

  }

# plot single day
plot.day.station2 <- function(jdf){
  
  # if(class(jdf$loc_time)[1]!='POSIXct'){
  #   jdf$loc_time <- as.POSIXct(jdf$loc_time)
  # }
  
  colors = c("Station" = "forestgreen", 
             "HRRR DSWRF" = "grey20",
             "HRRR DSWRF cosZ" = "grey60",
             "EKDS model cosZ" = "tomato3", 
             "Clearsky cosZ" = "deepskyblue1")
  
  shapes_ = c("Station" = 17, 
              "HRRR DSWRF" = 3,
              "HRRR DSWRF cosZ" = 1,
              "EKDS model cosZ" = 1, 
              "Clearsky cosZ" = 1)
  jdf %>% 
    # mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
    # mutate(station = visin * 0.67 + nirin * 0.33) %>%
    mutate(station = visin) %>%
    # mutate(station = visin + nirin) %>%
    ggplot(aes(x = loc_time)) + 
    # Line plot for 'Station' with correct color mapping
    geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
    
    # HRRR DSWRF with smoothing (ignoring NA's)
    geom_point(aes(y = dswrf, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.4, shape=3) +
    geom_line(data = subset(jdf, !is.na(dswrf)), aes(y = dswrf, colour = "HRRR DSWRF"), size = 0.75, alpha = 0.6, 
              linetype='dotted') +
    
    geom_point(aes(y = dswrf_z, colour = "HRRR DSWRF cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(dswrf_z)), aes(y = dswrf_z, colour = "HRRR DSWRF cosZ"), size = 0.5, alpha = 0.6) +
    
    # Point and smoothing layers for 'kn_mod', 'k_mod', and 'clearsky' with specific colors
    geom_point(aes(y = kn_mod_z, colour = "EKDS model cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(kn_mod_z)), aes(y = kn_mod_z, colour = "EKDS model cosZ"), size = 0.5, alpha = 0.6) +
    
    geom_point(aes(y = clearsky, colour = "Clearsky cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(clearsky)), aes(y = clearsky, colour = "Clearsky cosZ"), size = 0.5) +
    scale_color_manual(name = "", values = colors) +
    scale_x_datetime(labels = scales::date_format("%H:%M"), 
                     breaks = scales::date_breaks("5 hours"), timezone = "US/Mountain") +
    xlab('') + 
    ylab(expression(Irradiance ~ (W ~ m^-2))) + 
    theme_classic() + 
    theme(legend.position = c(0.9, 0.9),  # Position legend inside the plot area
          plot.margin = margin(10, 20, 30, 30))  # Adjust plot margins
  
}

# plot single day
plot.multiday.station <- function(jdf){
  
  if(any(class(jdf$loc_time)!='POSIXct')){
    jdf$loc_time <- as.POSIXct(jdf$loc_time)
  }
  
  colors = c("Station" = "forestgreen", 
             "HRRR DSWRF" = "grey20",
             "HRRR DSWRF cosZ" = "grey60",
             "EKDS model cosZ" = "tomato3", 
             "Clearsky cosZ" = "deepskyblue1")
  
  shapes_ = c("Station" = 17, 
              "HRRR DSWRF" = 3,
              "HRRR DSWRF cosZ" = 1,
              "EKDS model cosZ" = 1, 
              "Clearsky cosZ" = 1)
  
  jdf %>% 
    mutate(station = visin) %>%
    ggplot(aes(x = loc_time)) + 
    # Station line plot
    geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
    # HRRR DSWRF with smoothing (ignoring NA's)
    geom_point(aes(y = dswrf, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.4, shape=3) +
    geom_line(data = subset(jdf, !is.na(dswrf)), aes(y = dswrf, colour = "HRRR DSWRF"), size = 0.75, alpha = 0.6, 
              linetype='dotted') +
    # HRRR DSWRF with smoothing (ignoring NA's)
    geom_point(aes(y = dswrf_z, colour = "HRRR DSWRF cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(dswrf_z)), aes(y = dswrf_z, colour = "HRRR DSWRF cosZ"), size = 0.5, alpha = 0.6) +
    # EKDS model with smoothing (ignoring NA's)
    geom_point(aes(y = kn_mod_z, colour = "EKDS model cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(kn_mod_z)), aes(y = kn_mod_z, colour = "EKDS model cosZ"), size = 0.5, alpha = 0.6) +
    # Clearsky with smoothing (ignoring NA's)
    geom_point(aes(y = clearsky, colour = "Clearsky cosZ"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf, !is.na(clearsky)), aes(y = clearsky, colour = "Clearsky cosZ"), size = 0.5) +
    # Color and formatting
    scale_color_manual(name = "", values = colors) +
    scale_x_datetime(labels = scales::date_format("%b %d"),
                     breaks = scales::date_breaks("1 day")) +
    xlab('') + ylab(expression(Irradiance ~ (W ~ m^-2))) + 
    theme_classic() + 
    theme( axis.text.x = element_text(angle = 45, hjust = 1),
           # legend.position = c(0.9, 0.9),
          plot.margin = margin(10, 20, 30, 30))
  
}

plot.multiday.diff <- function(jdf){
  
  if(any(class(jdf$loc_time)!='POSIXct')){
    jdf$loc_time <- as.POSIXct(jdf$loc_time)
  }
  
  # plot diff
  colors = c( "Station" = "forestgreen", 
              "EKDS model" = "tomato3",
              "HRRR DSWRF" = "grey20",
              "Clearsky" = "deepskyblue1")
  jdf2 = jdf %>% 
    mutate(station = visin) %>%
    mutate(dswrf_bias = station - dswrf) %>% 
    mutate(kn_mod_bias = station - kn_mod_z) %>% 
    mutate(clear_bias = station - clearsky) 
  jdf2 %>% ggplot(aes(x = loc_time)) + 
    geom_line(aes(y = station, colour = "Station"), size = 1.5, alpha = 0.4) +
    # variables
    geom_hline(yintercept = 0, linetype='dashed', color = 'grey60',size = 0.5) +
    geom_point(aes(y = dswrf_bias, colour = "HRRR DSWRF"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf2, !is.na(dswrf_bias)), aes(y = dswrf_bias, colour = "HRRR DSWRF"), size = 0.5, alpha = 0.6) +
    geom_point(aes(y = kn_mod_bias, colour = "EKDS model"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf2, !is.na(kn_mod_bias)), aes(y = kn_mod_bias, colour = "EKDS model"), size = 0.5, alpha = 0.6) +
    geom_point(aes(y = clear_bias, colour = "Clearsky"), show.legend = TRUE, size = 2, alpha = 0.2) +
    geom_line(data = subset(jdf2, !is.na(clear_bias)), aes(y = clear_bias, colour = "Clearsky"), size = 0.5) +
    scale_color_manual(name = "", values = colors) +
    scale_x_datetime(labels = scales::date_format("%b %d"),
                     breaks = scales::date_breaks("1 day")) +
    xlab('') +     ylab(expression(Model ~ Bias ~ (W ~ m^-2))) + 
    theme_classic() + 
    theme( axis.text.x = element_text(angle = 45, hjust = 1),
           # legend.position = c(0.9, 0.9),
           plot.margin = margin(10, 20, 30, 30))
}

map.plotday.gif <- function(grib_date_path, plot_h_level = c(1,2,3), which.hrrr = 1){
  #
  # 1. plot h
  # 2. plot hsw
  # 3. plot hsw3
  grib_list = list.files(grib_date_path, pattern = 'wrfsfcf06.grib2', full.names = T)
  # bounds
  topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
  hsh = hillShade(terrain(dem, opt='slope',unit='radians'), terrain(dem, opt='aspect',unit='radians'), 
                  angle = 25, direction = 145)
  for (g in grib_list){
    #
    print(basename(g))
    # read in h
    read.hrrr(g, dpoly, which.hrrr = which.hrrr) # returns: h, dtime (UTC), jdate 
    if(plot_h_level < 2){
      # plot h
      plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE, main = dtime)
      plot(h[['dswrf']], add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
      plot(h[['dswrf']], legend.only=TRUE, col=heat.colors(100, 0.5),
           legend.width=2.5, legend.shrink=0.7,
           axis.args=list(cex.axis=1.2),
           legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
      legend("bottomleft", c("DSWRF (h)",
                             paste(" grib_tz:", strsplit(basename(g), "\\.")[[1]][2]),
                             format(dtime, tz="MST",usetz=TRUE)), bty='n')

    }else if(plot_h_level <3){
      sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
      # plot hsw
      plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE, main = dtime)
      plot(hsw, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
      plot(hsw, legend.only=TRUE, col=heat.colors(100, 0.5),
           legend.width=2.5, legend.shrink=0.7,
           axis.args=list(cex.axis=1.2),
           legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
      legend("bottomleft", c("DSWRF (hsw)",
                             paste(" grib_tz:", strsplit(basename(g), "\\.")[[1]][2]),
                             format(dtime, tz="MST",usetz=TRUE)), bty='n')
    }else{
      clear.sky.sw(dem, dpoly, dtime, savepath=NULL) # returns sw3
      hrrr.sw.topo(hsw, cld, dem, savepath=NULL) # hsw3
      # plot hsw0-3
      plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE, main = dtime)
      plot(hsw3, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
      plot(hsw3, legend.only=TRUE, col=heat.colors(100, 0.5),
           legend.width=2.5, legend.shrink=0.7,
           axis.args=list(cex.axis=1.2),
           legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
      legend("bottomleft", c("DSWRF (hsw3)",
                             paste(" grib_tz:", strsplit(basename(g), "\\.")[[1]][2]),
                             format(dtime, tz="MST",usetz=TRUE)), bty='n')
      
    }
    
    rm(list=setdiff(ls(), c('dem', 'viewf', 'dpoly', 'hsh', 'grib_list', 'plot_h_level', 'which.hrrr') ) )
  
  }
  # rm(list=lsf.str())
}



# END PLOT FUNCTIONS
###























# # # # # # # # # # # # # # # # # # # # # # # #
# GENERAL - helper functions


# # to do
# k_modify_clouds <- function(k, Gh, cld_lyr){
#   # Modify k for cloud layer
#   Dh = 0 + 1*Gh
#   return(Dh)
# }

###############################################
# generate k values from HRRR
sw.splitter <- function(hsw, dem, sw0, Id0, zenith, jdate){
  #
  # split dswrf into Ib and Id based on diffuse fraction (k)
  # input:
  #     hsw - elevationally distributed HRRR dswrf
  #     sw0 - Ghc - modeled clearsky global horizontal irradiance
  #     Id0 - Dhc - modeled clearsky diffuse irradiance
  #     zenith - solar zenith angle
  #     jdate - julian date
  #     
  # returns:
  # k - Arias-Ruiz 2010 Boulder G2 
  k <<-  diffusef_Lk2_boulder(IG = hsw, zenith, I0 = exo_irr(jdate), m_elv = dem)
  #
  # newk - Ariz-Ruiz and Gueymard 2024
  new_k <<-  diffusef_gisplit4_cloudless(Gh = hsw, 
                                         Ghc = sw0, Dhc = Id0, E0h = exo_irr(jdate)) 
}

###############################################
diffusef_gisplit4_cloudless <- function(Gh, Ghc, Dhc, E0h) {
  # Function to calculate K using the given equation
  # After Ruiz-Arias and Gueymard, 2024 (GISPLIT)
  # inputs:
  # Gh:  GHI or global horizontal irradiance 
  # Ghc: the modeled GHI under an ideal cloudless sky,  
  # Dhc: the modeled DIF under the same ideal cloudless sky 
  #
  KT = Gh/E0h   # KT = Gh/E0h, where E0h is the extraterrestrial horizontal solar irradiance
  Kcs = Gh/Ghc  # Kcs = Gh/Ghc - clear sky index
  Kds = Dhc/Ghc # Kds = Dhc/Ghc - clear sky diffuse fraction
  # Parameters for airmass -  f0; f1; f2; f3; f4
  m <- airmass(zenith) # simple airmass for now
  fpars <- params_am_cloudless(m)
  # argument inside the exponential
  exponent <- as.numeric(fpars['f0']) + as.numeric(fpars['f1']) * KT + as.numeric(fpars['f2']) * KT + as.numeric(fpars['f3']) * Kcs + as.numeric(fpars['f4']) * Kds
  # K using the logistic function
  K <- 1 / (1 + exp(exponent))
  return(K)
}

###############################################
params_am_cloudless <- function(m){
  # tables in Ariaz-Ruiz and Gueymar 2024
  # https://pdf.sciencedirectassets.com/271459/1-s2.0-S0038092X24X0002X/1-s2.0-S0038092X24000574/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjENL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDlRw39ol1Pg5Ci%2FX5mIS2h%2BkV4%2BTTRezEEy%2BO8kqKYowIgIvAWA%2FJH8pW9o%2FnBwl5VSvrmBYlPt1LAYZZ6HGuMNugqswUIexAFGgwwNTkwMDM1NDY4NjUiDGAk8oYvCAdhtgemRSqQBZXh2ERfxbHH003wnxQx42e%2FfMW9Mp0vZ0%2FBitcOCIwkW20AxGUePDQxpmBj%2FAFMyvc6QfH1Iu5e6IXHUmZVFB4dotF0J5MvbMALV682am871zZ%2FTN18LQeLjfhcI%2FVujH9l5PmwM6eElff0q4UAWpr6cyadFJEI0idQSRGvpC8KewvOsg8H59mYbp44e7INXolX17y5jF83pM4IkBJ6ALbIrzwSaj%2B2QN5KPAfsjpR3Hs0OLf60Uqdwvu7zey4f35alJABgF9JEtSF%2FOEoJw3hb6VdYZz3HqkBVovgr9nfDDAIJ4QbNGhLqKwWcHB8tnKZ2m%2BBKkMJgxzNzn1apXuwjXpz8ExAd%2FFhDkN51k3GVNMzw85xiKdXxjyet3QF0VUYSm8WEg75Ft%2F9PCPxu%2FQ%2Fcy%2Bu4DlKZB2rc0L1hub13JBqvVe6H28RBUGKZeBGQKOZ87Czi0xgMl1BHnT4qcPP4yitnFOCnMmLm5w9VwzjV0NFYuRS%2BI3PGoyU%2BDXTghqgSCyyOACwCjQDMUJ%2FgA5Z%2BvV2%2FCUFwD01arcCAiSRCHgDGLRES9z3kcxLvQp7MfkMNc5eZnLhAhdQKdlSReHnOlKC7JPIVzEijwZYS5aIpsOatFzxo6CjcrHUE0LqAF5VOVFvhuVqJpjCqgyMygM3ZV5WsMRMu%2FfUfr%2FZ8BZlDZE86RQPq0xN%2BlSvbpo96wBQ%2B1msQKDnVSqALqiQSSOy6dHc9OAcqUX%2B8C55r%2F%2B9MqjhUDC7E8jo4T1l3LTv8U7EFMHpWkEgitmCycsQPPRjwkJb9k9s2bXeEfTnTH%2FCOgLAQHqDndlL6IGVMakiaRsO0M8wmg3t%2FWkeF6HwFx2WHayXhzdprSk%2FHHFOiP5mMMMaCqLoGOrEB9RRhQ2b2uXkEpRINHtzBohzafvIhgU63uZ9OkfqyV8NIwAZ2Ao5eFLCrwGqhxtki66wt45HfQflEH4VrOwhjwDtXjAd4%2Bq6ZnFD%2FMYy4n3m1Fgd5%2B4hvRan4d5yTWOsWg6zL1XpTXyZh75eUdz9zut87wnpsJaIVKaZiVD2IIo1l7su49lmVHtYbBbn192xkMULoLAO0G3VJn9mMhV0WGZaoZv3M%2FszzGzdLykdnj5yO&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241129T183348Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY7KZCLWXM%2F20241129%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e7a306d45318c8aecafec449f8f97edfde8f6c7eeb42a521d24b2828297a2f4d&hash=02eab031a706025fc9af0606f4051a8c904d434c78187eb7a7137b1c6d54ec27&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0038092X24000574&tid=spdf-e4c1d2bc-0729-47e2-b327-3d46784d55d7&sid=dc8d14b12521c24078197d9422b802d85440gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f155d0a02530d5e54095d&rr=8ea4b02f8ab327d4&cc=us
  am_vals1 = c(1.0, 1.2115276586285884, 1.4677992676220695, 1.7782794100389228, 2.154434690031884, 2.610157215682537,
               3.1622776601683795, 3.831186849557287, 4.641588833612778, 5.623413251903491, 6.812920690579611, 8.254041852680183)
  m_match = data.frame( range = 1:12, am_lower = am_vals1,
                        am_upper = c(am_vals1[2:12], 10.0) )
  # params in 'gisplit_for_obs_01min_caelus_reg'
  table <- list(
    c(-1.1927111865238231, 1.0626052976611577, 3.192159399906436, 0.5800276621224962, -4.287378904381701),
    c(-1.4639374536744845, 2.584489765207779, 1.6100156006572253, 0.8098018976716924, -3.8008175914191624),
    c(-1.0368390353097843, 4.340858468717979, -1.8258808978843784, 1.6161591699273021, -3.7528424477351563),
    c(-0.623278706264383, 1.4685241946998655, 0.09355086127214213, 1.8495451023720764, -3.709492430104688),
    c(-0.46382678817272655, 1.4320638191469477, -0.42303090290827317, 2.00610060011075, -3.625812045518225),
    c(-0.378255971190094, -0.11293277352589735, 0.8487339289470736, 2.014385345212601, -3.462603304681434),
    c(-0.33783820336303555, -0.7011053769485842, 1.2365297466014626, 2.000446054790625, -3.320745470520213),
    c(-0.4699887060267603, 0.7284885905251159, -0.5518540448948532, 2.3285365753267797, -3.445325299775386),
    c(-0.8324966659219617, -3.323309518614715, 3.086877934795852, 2.9274838836621218, -3.7210252379651516),
    c(-1.1212458179161704, -3.7681245668951306, 3.050348260687978, 3.54244975958997, -4.179148820782683),
    c(-1.4971606261115933, -4.554618596896216, 3.4873825693735983, 4.164518683801098, -4.666295680462231),
    c(-1.5034778218634748, -3.5797469138080786, 1.879192604435627, 4.692566360761011, -5.44832996369178)
  )
  # convert to  matrix
  table_matrix <- do.call(rbind, table)
  # convert to dataframe
  table_df <- as.data.frame(table_matrix)
  colnames(table_df) <- c("f0", "f1", "f2", "f3", "f4")
  table_df$am_range <- 1:12
  # extract correct parameters based on airmass
  m2 = m_match$range[m_match$am_lower <= m & m_match$am_upper > m]
  return(table_df[m2,1:5])
}

###############################################
localG2 <- function(kt, m) {
  # Local G2 model (Boulder, CO) - from Arias-Ruiz 2010
  # params from text: https://pdf.sciencedirectassets.com/271098/1-s2.0-S0196890410X00036/1-s2.0-S0196890409004695/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEHwaCXVzLWVhc3QtMSJIMEYCIQC6Azl0wytBnGRHsMSYo8eDLC4UnkQsdKLRSUwpuyYDvgIhALsVStnviro5FrOeXkJnmuqo7sMZBs4%2BOI%2BUaSOIRhHeKrIFCCUQBRoMMDU5MDAzNTQ2ODY1IgxioNKhpG%2BL%2FDBFCisqjwXVClh0O%2B0KfcUoC%2FuasZdybQnzuxIH94YzUhZ90r8kJi%2FzCEWT9gAQorVPBnYRo1OxHP%2FAa5F%2BZp8znPglppziRLAUxPIx%2BbqNxZkasb%2B%2Bv6KCdHzBDjSEj%2FzBFdqTCQ0jfXblFmEcIwL21cgtR4DZzRCgApJQpotxbvCvirXQDIfUyWREWUQBuM3nr1lETVlau4jvcNrFxlzyTM5OJKWYE8tKJoTF5IeWa2biX2LeHIzMz8dd65CapyysRZvJC2f8Qvitt2KA0a%2FhcWAPZ7nesyp6LLeYr2fizJdDrIYWi7lVRmz2Pv4KwrSTQRK48CUaaf7BpJwNN4Pj%2BYNWL0L5vKRLvVjR3tiJw5LI1M1yqKPY3JK5quig7s2eGtXey%2BDwkKT9vwUufXj%2FgoV5VSWiJGJFKgMriTosqsskyzojWzChKEqJ5e7faoEt5oS1Nxjy59BDKG7bOQmXGXdiX9FgTO0caYUk2225fyVrrR8jqRFj7mfEvXbmnSqqgBJFyrAxcHr2gotDO7dBTEMVLChIxZDxIrsMUAK6aieCeaOEUATLWT6zzNUcoznmlR3ahY16Z0ac5tKoyuNcGQXM7tqW9nQ%2BEfzSJGVHyzJaCVVZ019WzxgMX0mPb184IaLFOe6Ngse03ui4ODou6JXuoecOVbKCly8pTW9kWW5D0hsiU0NJf4EOJbp8%2Bt0H1Fu5kgEUa9oosZD%2FZ%2FHyrv4VepiLFWUlFZlI3cfUdK5fo52nZjoLvQ96rOs%2BwC4yg3qyBCNUVgcVNdK3paFvDCUx68T05OiqHQmArDGZDPb7imO98iP7HYy2Bf8ACVNKXm0B2YV0PltgwjuqxR3fSydBKGh8R78%2BYlFBcb%2FG9GtP6Dq9MIiElboGOrABM0H7OSU6aagORuIjQTt31Uegr%2F1G%2B3ncRUjf50ieMJ3yIszL7%2BkeRufJ3x%2FIRZILjG02d3yK7q2DuOUfpIMLj6kzzrJoRrWwsFTmoGbCbrCld%2FmF0JS7u4wOj7RhL3aa6Yz4izHgqxkU0gFkfQWrOgzxP0YnTeggPFMVXvs8HrcI%2FI8u%2FnAUMniX2a4RCxlZkSIYSSbVE2cVy8ETK08pHMLaE%2FciRr%2FUaGpeYJxDnCU%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241126T040419Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYZTSKULMN%2F20241126%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0d47411edac7df87f4b54d81da4665f3a317524af6db7eba11dcf6af4f9729a3&hash=6de1c4facec3d285a69f64febee99618e1bb28b9c69ec85e7261306f700654ee&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0196890409004695&tid=spdf-efeda5c7-5e71-44f6-970d-0408f2ae8751&sid=4afaa25e92dee841b98abc64158234a77a65gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f155d0a5b51090b505a00&rr=8e86fe65ead26a05&cc=us
  a0= 0.956; a1= 1.268; a2= 3.202; a3= -6.712; a4= 2.228; a5= -0.213; a6= 0.021; 
  # inner exponential expression
  inner_exp <- -exp(a2 + a3 * kt + a4 * kt^2 + a5 * m + a6 * m^2)
  # full function value
  k_value <- a0 - a1 * exp(inner_exp)
  return(k_value)
}

airmass <- function(zenith){
  # m - airmass simple
  m = 1/cos(zenith * pi / 180)
  return(as.numeric(m) )
} 

###############################################
elevational_airmass <- function(zenith, elevation) {
  # elevation-adjusted air mass based on Young (1987)
  # base air mass calculation
  base_m <- 1 / cos(zenith * pi / 180)  # Convert zenith to radians
  # elevation adjustment using the barometric formula for pressure ratio
  pressure_ratio <- (1 - (2.25577e-5 * elevation))^5.2559
  # adjusted air mass
  adjusted_m <- base_m * pressure_ratio
  return(adjusted_m)
}

###############################################
diffusef_Lk2_boulder <- function(IG, zenith_hr, I0 = 1367, m_elv = NULL) {
  # calculate diffuse fraction (k0) using clearness, airmass, and diffuse fraction models
  # clearness index (kt)
  # kt <- clearness(IG, zenith_hr, I0)
  kt = IG/I0
  # calculate airmass (m) 
  if (!is.null(m_elv)) {
    # airmass calculation if elevation is provided
    m <- elevational_airmass(zenith_hr, m_elv)
  } else {
    # airmass calculation if no elevation is provided
    m <- airmass(zenith_hr)
  }
  # diffuse fraction (k0) using clearness and airmass
  k0 <- localG2(kt, m)
  return(k0)
}

###############################################
exo_irr <- function(day_of_year) {
  # Solar constant (W/m^2)
  I_sc <- 1361
  # Calculate extraterrestrial solar irradiance using the formula
  I_0 <- I_sc * (1 + 0.034 * cos((2 * pi * day_of_year) / 365))
  return(I_0)
}

###############################################
insol.params <- function(use_hrrr=FALSE){
  visibility <<- 25 
  O3 <<- 340/1000 # DU = 0.01 mm (CO 325-340 DU)
  # https://disc.gsfc.nasa.gov/datasets?page=1&source=AURA%20OMI
  alphag <<- 0.35
  if(use_hrrr){
    # print("! currently unable to solve for rh and temp at sub grid-level")
    # stop()
    # RH <- array(h[['rh']])
    RH <<- array(resample(h[['rh']], dem))
    # tempK <- array(h[['tempC']]) + 273.15 #278.15 # 5 C or 41 F
    tempK <<- array(resample(h[['tempC']], dem)) + 273.15
  }else{
    RH <<- 40 # 60
    tempK <<- 273.15 #278.15 # 5 C or 41 F
  }
}

###############################################
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}

###############################################
slope.aspect <- function(dem, units = 'radians', neighbor = 8){
  s <- terrain(dem, opt='slope',unit=units,neighbors=neighbor)
  a <- terrain(dem, opt='aspect',unit=units,neighbors=neighbor)
  stk <- stack(s,a)
  return(stk)
}

###############################################
cos.slope <- function(zenith, azimuth_eq, aspect, slope, as_mat=FALSE){
  # returns a matrix of the cosine of the incident angle at a given moment
  exposures  = aspect - radians(180)
  cos_inc = acos((cos(radians(zenith)) * cos(slope)) +
                   (sin(radians(zenith)) * sin(slope) * cos(radians(azimuth_eq) - exposures)))
  if(as_mat){cos_inc = as.matrix(cos_inc)}
  # get rid of self shading values
  cos_inc[cos_inc > radians(90)] = radians(90)
  cos_inc = cos(cos_inc)
  return(cos_inc)
}

###############################################
view.factor <- function(dem_mat, dem, dem_res, elv_interval = 5, az_interval = 15){
  print("____Generating sky view factor_____")
  ptm <- proc.time()
  ELV = rev(seq(0, 90, elv_interval))
  AZI = seq(0, 345, az_interval)
  AZ = matrix(0,nrow=dim(dem_mat)[1]*dim(dem_mat)[2],ncol=length(AZI))
  for (vv in 1:length(AZI)){
    print(paste(vv, "of", length(AZI)))
    Z1 = matrix(0,nrow=dim(dem_mat)[1]*dim(dem_mat)[2],ncol=length(ELV))
    for (mm in 1:length(ELV)){
      sv = normalvector(ELV[mm],AZI[vv])
      sh <- doshade(dem_mat, sv, dl=dem_res)
      Z1[,mm] = as.array(sh)
    }
    AZ[,vv] = rowSums(Z1)/length(ELV)
  }
  VF_mat = matrix(rowMeans(AZ), nrow=nrow(dem), ncol=ncol(dem))
  #VF_dem <- make.raster(VF_mat, dem)
  print(proc.time() - ptm)
  return(VF_mat)
}