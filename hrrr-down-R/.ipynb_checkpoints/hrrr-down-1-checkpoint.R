#
# HRRR Downward Shortwave Radiation Elevational Downscaling
# R 4.03 Geospatial Packages (CHPC)
# HRRR v04
#
# # # # # # # # #
rm(list = ls());gc()

require(insol);require(raster);require(rgdal);require(rgeos)
require(ncdf4); require(sp)

# # # # # PATHS
# TOPO.nc FILE
# dem, mask, veg_type, veg_height, veg_k, veg_tau, Alec2 mask, sky_view_factor, terrain_config_factor, slope
dem = raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/topo.nc", varname = 'dem')
crs(dem) = "EPSG:32613"
plot(dem)
# sky view
sv = raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/topo.nc", varname = 'sky_view_factor')
plot(sv, col = grey(1:100/100))
slp = raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/topo.nc", varname = 'slope')
plot(slp, col = grey(1:100/100))
# DEM area
dpoly = as(extent(dem), 'SpatialPolygons')
crs(dpoly) = crs(dem)
# save for KML
# d.df <- data.frame( ID=1:length(dpoly), row.names = "1") 
# p <- SpatialPolygonsDataFrame(dpoly, d.df)
# pdpoly <- spTransform(p, "EPSG:4326")
# pdpoly@data$name <- "topo_nc_extent"
# write_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw'
# # write_path = '/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/'
# writeOGR(obj = pdpoly, dsn = file.path(write_path, "topo_boundary.kml"), layer="topo_boundary", driver="KML")
# list.files(write_path)

# # # # #
# SNOTEL
'/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/Irwin_radiation_data-2019-2022.xlsx'
idf = read.csv('/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/wy2022/Irwin/2022-Irwin.csv', row.names = NULL)
head(idf) # max 1300 min 0 - no out when no in 
names(idf) <- c('time', 'temp', 'visin', 'nirin', 'visout', 'nirout')

idf$visin[idf$visin > 1300] = 1300; idf$visin[idf$visin < 0] = 0
idf$nirin[idf$nirin > 1300] = 1300; idf$nirin[idf$nirin < 0] = 0

as.Posixct(idf$time[1])
as.POSIXct(idf$time[1], "%m/%d/%Y %H:%M", tz = "UTC")
idf$time = as.POSIXct(idf$time, "%m/%d/%Y %H:%M", tz = "UTC")

idf %>% mutate(dswrf = visin + nirin) %>% 
  # mutate(time = as.POSIXct(idf$time[1], "%m/%d/%Y %H:%M", tz = "UTC")) %>% 
  ggplot(aes(x = time, y = dswrf)) + geom_line() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_x_datetime(labels = date_format("%Y-%m-%d")) +
  xlab("")

# Step 1: Define the start and end dates for the range
start_date <- as.POSIXct("2022-03-26 00:00:00") # 19 for two weeks
end_date <- as.POSIXct("2022-04-01 23:59:59")

# Step 2: Subset the data based on the date range
idf2 <- idf[idf$time >= start_date & idf$time <= end_date, ]
  
idf2 %>% mutate(dswrf = visin + nirin) %>% 
  ggplot(aes(x = time, y = dswrf)) + geom_line() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_x_datetime(labels = scales::date_format("%Y-%m-%d")) +
  xlab("")




# new 30 m DEM
list.files('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw/')
save_base_name = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw'
ls_r <- list.files('/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/toponc_ERW/srtm_tiles/',
                   pattern = ".tif", full.names = T)
rl <- lapply(ls_r, raster)
cat("Reprojecting and saving SRTM DEM mosaic...")
# ddproj = projectRaster(dem, crs = crs("EPSG:4326"))
pdpoly <- spTransform(dpoly, "EPSG:4326") 
shape <- buffer(pdpoly, 0.02) # lat long buffer
# dem3 = do.call(merge, c(rl, tolerance = 1)) %>% crop(.,shape) # for testing buffer
# save
stop()
do.call(merge, c(rl, tolerance = 1)) %>% 
  crop(.,shape) %>% 
  projectRaster(., crs = crs(dem), filename = file.path(save_base_name,"erw_srtm30_merged.grd"), overwrite=T)

# need 64 GB memory
# do.call(merge, c(rl, tolerance = 1)) %>%
#   projectRaster(., crs = crs(dem)) %>%
#   crop(.,spTransform(shape, proj4string(dpoly) ), filename = file.path(save_base_name,"erw_srtm30.grd"), overwrite=T)

# dpath2 = "/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw/erw_srtm30_merged.grd"
dpath2 = "/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/erw/erw_srtm30.grd"
file.exists(dpath2)
dem2 = raster(dpath2)
plot(dem2)

res(dem2)
viewf = view.factor(as.matrix(dem2), dem2, 90) # 24.1
viewf2 = make.raster(viewf, dem2)
plot(viewf2, col = grey(1:100/100))
# writeRaster(viewf2, file.path(save_base_name,"erw_sv24.grd"))

vv = raster(file.path(save_base_name,"erw_sv24.grd"))

plot(shape)
plot(pdpoly,add=T)
plot(dem2, add=T)

s_a = slope.aspect(dem2)
sdv = stack(dem2, s_a[[1]], s_a[[2]], viewf2)
names(sdv) = c('dem', 'slp', 'aspect', 'sky_view_factor')
# writeRaster(sdv, file.path(save_base_name,"erw_topo.grd"), overwrite=T)
# writeRaster(sdv3, file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/grand-mesa/dem/",
#                            "gm_topo2.grd"), overwrite=T)


"+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs"

sdv2 = brick("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/grand-mesa/dem/gm_topo.grd")
sdv3 = projectRaster(sdv2, crs="+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs") 
sdv3 = crop(sdv3, ppp)

r = raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/grand-mesa/dem/n39_w109_3arc_v2.tif")




# ATWATER
dest_file = "/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/data_GSLB"
dpath_jord = '/uufs/chpc.utah.edu/common/home/skiles-group2/otto/GSLB_iSnobal_data/Outputs/Jordan/topo.nc'
# dem, mask, veg_height, veg_k, veg_tau, veg_type, sky_view_factor, terrain_config_factor, slope
demj = raster(dpath_jord, varname = "dem") # huge
viewfj = raster(dpath_jord, varname = "sky_view_factor") # huge
plot(demj)
plot(viewfj)


# nc = nc_open(file.path(save_base_name,"erw_topo.nc"), write = TRUE)
# names(nc$var)
# ncdf4::ncvar_put(nc, 'terrain', c('dem', 'slp', 'aspect', 'sky_view_factor'))

# shape
irwin_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/site-locations/irwin_snotel.json'
('/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/Snotel/site-locations/snotel_sites_3x3.json')
con = file(irwin_path)
line = readLines(con)
lon0 = trimws(grep("lon", line, value = T))
lat0 = trimws(grep("lat", line, value = T))
x = as.numeric(strsplit(strsplit(lon0, ": ")[[1]][2], ',')[[1]][1])
y = as.numeric(strsplit(strsplit(lat0, ": ")[[1]][2], ',')[[1]][1])
# points
pts = cbind(x, y)
pts <- SpatialPoints(coords = pts, proj4string = CRS(proj4string(dem)) )
plot(pts, add=T,col='red', cex=1, pch=20)
pb = buffer(pts, 3.5e3) # 3500 m surrounding cells
plot(extent(pb), add=T, col='red')


# # # # #
# HRRR DATA
# smrf_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/erw_isnobal/wy2022/erw_hrrr_solar/run20220325/net_solar.nc'
smrf_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/erw_isnobal/wy2022/erw_hrrr_solar/run20220401/net_solar.nc'
net_solar = stack(smrf_path, varname="net_solar")
plot(net_solar[[16]])
# DSWRF, illumination_angle, albedo_ir, albedo_vis, albedo, net_solar
# cdo -f nc copy file.grb file.nc
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/hrrr.t13z.wrfsfcf06.grib2'

# nc = nc_open(grib_path)
# vars
# Geopotential height [gpm]
# Temperature [C]
# Relative humidity [%]
# u-component of wind [m/s]
# v-component of wind [m/s]
# 06 hr Total precipitation [kg/(m^2)]
# 01 hr Total precipitation [kg/(m^2)]
# Total cloud cover [%]
# Total cloud cover [%]
# Downward short-wave radiation flux [W/(m^2)]
# Visible Beam Downward Solar Flux [W/(m^2)]
# Visible Diffuse Downward Solar Flux [W/(m^2)]

# HRRR
hrrr = stack(grib_path)
hrrr
names(hrrr) <- c('gmp', 'tempC', 'rh', 'u-wind', 'v-wind', 'prec06hr','prec01hr', 'cld0','cld', 'dswrf','vbdsf', 'vddsf')
plot(hrrr[[10]])

plot(hrrr[['dswrf']]/hrrr[['vbdsf']])




# clip and crop
dpoly_proj = spTransform(dpoly, proj4string(hrrr) )
plot(dpoly_proj, add=T)
# reprojected 
hrc = crop(hrrr, dpoly_proj) #, snap='out')
# hrc = crop(hrrr, buffer(dpoly_proj, 3.5e3)) # + 1 grid cell
# h = projectRaster(hrc, dem) # change from 3000 m to 50 m
h = projectRaster(hrc, crs = crs(dem) )
# solar
plot(h[[10]], col = heat.colors(100))
plot(dpoly, add=T)
plot(pts, add=T,col='blue', cex=1, pch=20)

plot(values(h[[10]])~values(h[[1]])) # more of exponential?
extract(h, pts)

# zoom
hh = crop(h, pb)
dim(hh)
plot(hh[[10]])
plot(pts, add=T,col='red', cex=1, pch=20)

d = crop(dem, hh)
dim(dem)
plot(d)
plot(pts, add=T,col='red', cex=1, pch=20)

# single grid cell
ir_hgrid = crop(hh, buffer(pts, 1.5e3))[[10]]
dgrid = crop(d, ir_hgrid)
plot(ir_hgrid)
plot(dgrid, add=T) 
plot(pts, add=T,col='red', cex=1, pch=20)

# not possible
# second_derivative = 2 * (G_z - Gz0) / var_z


# ELV DOWNSCALE
hd <- raster(dem)
hd[] <- 1
hdd  = disaggregate(h[['dswrf']], fact = 60) # 300/50
hdd = resample(hdd, hd)
hdd2 = hdd

# same resolution - now use slope to correct
plot(hdd)

# relationship
plot(values(h[['dswrf']])~values(h[['gmp']]) )
dat2 = data.frame(x = values(h[['gmp']]), y =values(h[['dswrf']])) %>% na.omit()
model <- lm(y ~ x, dat2)
intercept <- coef(model)[1]  # Intercept (b0)
slope <- coef(model)[2]     # Slope (b1)
# predicted_y <- intercept + slope * new_x
hrd1 <- intercept + slope * dem
hrnn = projectRaster(hrc[['dswrf']], dem) # change from 3000 m to 50 m
hrbil = projectRaster(hrc[['dswrf']], dem, method='bilinear') # or ngb
plot(hrd1)
plot(hrnn)
plot(hrbil)
extract(hrd1, pts)
extract(hrnn, pts)
extract(hrbil, pts)
# show for zoomed in gridcell up close

# other things
# improve elevational distribution of solar
# get new DEM 30X30 expanded
# export topo for zoomed in focus on IRWIN
# better Kt and K separation


# estimate K 
exo_irr <- function(day_of_year) {
  # Solar constant (W/m^2)
  I_sc <- 1361
  
  # Calculate extraterrestrial solar irradiance using the formula
  I_0 <- I_sc * (1 + 0.034 * cos((2 * pi * day_of_year) / 365))
  
  return(I_0)
}
grib_hr <- gsub(".*t([0-9]+)z.*", '\\1', basename(grib_path) )
day_date <-  as.Date(strsplit(basename(dirname(grib_path)), "\\.")[[1]][2], "%Y%m%d")
jdate <- as.numeric(format(day_date, "%j"))

kt = hrrr[['dswrf']]/exo_irr(jdate)
plot(kt)
plot(dpoly_proj,add=T)



 

# INSOLATION MODEL
# # # # # # # # # # #
# timestamp
grib_hr <- gsub(".*t([0-9]+)z.*", '\\1', basename(grib_path) )
day_date <-  as.Date(strsplit(basename(dirname(grib_path)), "\\.")[[1]][2], "%Y%m%d")
moment_date_char <- paste0(as.character(day_date), " ", grib_hr, ":00:00")
dtime <- as.POSIXct(moment_date_char, tz = "UTC") # UTC +1 for grid

# SOLAR
# sun position and vector
jd = JD(dtime)

# pdpoly <- spTransform(dpoly, "EPSG:4326")
# lat_lon <- rowMeans(pdpoly@bbox)[2:1]
demh = h[[1]] 
ddproj = projectRaster(demh, crs = crs("EPSG:4326"))
lat_lon <- c(round((ddproj@extent@ymax + ddproj@extent@ymin)/2,5), round((ddproj@extent@xmax + ddproj@extent@xmin)/2,5))


sv = sunvector(jd,lat_lon[1],lat_lon[2],-7) 
sp1=sunpos(sv)
# daylight hours (zenith <= 90)
sp=sp1[which(sp1[,2]<=90),]
sv=sv[which(sp1[,2]<=90),]
# solar noon
zenith = sp[2]
azim = sp[1]
azimuth_eq = azim -180

#### HELPER FUNCTIONS
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}
slope.aspect <- function(dem, units = 'radians', neighbor = 8){
  s <- terrain(dem, opt='slope',unit=units,neighbors=neighbor)
  a <- terrain(dem, opt='aspect',unit=units,neighbors=neighbor)
  stk <- stack(s,a)
  return(stk)
}
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


# rerun shade
# make.shade(zenith, azim)
# if(file.exists(d_mat)){d_mat <<- as.matrix(demh)}
d_mat <- as.matrix(demh)
dl = res(h)[1]
sh0 <- doshade(d_mat, sv, dl=dl) # does not make sense for course grid
sh0 <- make.raster(sh0, demh) 
# zenith and incident angles
s_a <- slope.aspect(demh)
cos_inc <- cos.slope(zenith, azimuth_eq, aspect = s_a[[2]], slope = s_a[[1]])

# SOLAR NOON INSOLATION 
visibility <- 25 
RH <- 40 # 60
tempK <- 273.15 #278.15 # 5 C or 41 F
O3 <- 340/1000 # DU = 0.01 mm (CO 325-340 DU)
alphag <- 0.35
height <- array(demh)
# RH <- array(h[['rh']])
# tempK <- array(h[['tempC']]) + 273.15 #278.15 # 5 C or 41 F
Idirdif = insolation(zenith,jd,height,visibility,RH,tempK,O3,alphag)
Ib = matrix(Idirdif[,1],nrow=nrow(demh),ncol=ncol(demh), byrow = T)
Id = matrix(Idirdif[,2],nrow=nrow(demh),ncol=ncol(demh), byrow=T)
# MAKE RASTER
# cosi = make.raster(cos_inc0,d)
Ib0 = make.raster(Ib, demh)
Id0 = make.raster(Id, demh) 

# sw1 = Ib0 * cos_inc #* sh0 # direct beam only
sw1 = Ib0 * cos(radians(zenith))
sw2 = sw1 + Id0 # direct and diffuse
sw3 = Ib0 + Id0
# sw20 = sw2

plot(sw2)
plot(stack(sw3,h[['dswrf']]))
plot(sw3/h[['dswrf']])

plot(values(sw2)~values(h[['gmp']]))
plot(values(h[['dswrf']])~values(h[['gmp']]))
#
#








# Local G2 model (Boulder, CO) - from Arias-Ruiz 2010
localG2 <- function(kt, m) {
  # params from text
  # https://pdf.sciencedirectassets.com/271098/1-s2.0-S0196890410X00036/1-s2.0-S0196890409004695/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEHwaCXVzLWVhc3QtMSJIMEYCIQC6Azl0wytBnGRHsMSYo8eDLC4UnkQsdKLRSUwpuyYDvgIhALsVStnviro5FrOeXkJnmuqo7sMZBs4%2BOI%2BUaSOIRhHeKrIFCCUQBRoMMDU5MDAzNTQ2ODY1IgxioNKhpG%2BL%2FDBFCisqjwXVClh0O%2B0KfcUoC%2FuasZdybQnzuxIH94YzUhZ90r8kJi%2FzCEWT9gAQorVPBnYRo1OxHP%2FAa5F%2BZp8znPglppziRLAUxPIx%2BbqNxZkasb%2B%2Bv6KCdHzBDjSEj%2FzBFdqTCQ0jfXblFmEcIwL21cgtR4DZzRCgApJQpotxbvCvirXQDIfUyWREWUQBuM3nr1lETVlau4jvcNrFxlzyTM5OJKWYE8tKJoTF5IeWa2biX2LeHIzMz8dd65CapyysRZvJC2f8Qvitt2KA0a%2FhcWAPZ7nesyp6LLeYr2fizJdDrIYWi7lVRmz2Pv4KwrSTQRK48CUaaf7BpJwNN4Pj%2BYNWL0L5vKRLvVjR3tiJw5LI1M1yqKPY3JK5quig7s2eGtXey%2BDwkKT9vwUufXj%2FgoV5VSWiJGJFKgMriTosqsskyzojWzChKEqJ5e7faoEt5oS1Nxjy59BDKG7bOQmXGXdiX9FgTO0caYUk2225fyVrrR8jqRFj7mfEvXbmnSqqgBJFyrAxcHr2gotDO7dBTEMVLChIxZDxIrsMUAK6aieCeaOEUATLWT6zzNUcoznmlR3ahY16Z0ac5tKoyuNcGQXM7tqW9nQ%2BEfzSJGVHyzJaCVVZ019WzxgMX0mPb184IaLFOe6Ngse03ui4ODou6JXuoecOVbKCly8pTW9kWW5D0hsiU0NJf4EOJbp8%2Bt0H1Fu5kgEUa9oosZD%2FZ%2FHyrv4VepiLFWUlFZlI3cfUdK5fo52nZjoLvQ96rOs%2BwC4yg3qyBCNUVgcVNdK3paFvDCUx68T05OiqHQmArDGZDPb7imO98iP7HYy2Bf8ACVNKXm0B2YV0PltgwjuqxR3fSydBKGh8R78%2BYlFBcb%2FG9GtP6Dq9MIiElboGOrABM0H7OSU6aagORuIjQTt31Uegr%2F1G%2B3ncRUjf50ieMJ3yIszL7%2BkeRufJ3x%2FIRZILjG02d3yK7q2DuOUfpIMLj6kzzrJoRrWwsFTmoGbCbrCld%2FmF0JS7u4wOj7RhL3aa6Yz4izHgqxkU0gFkfQWrOgzxP0YnTeggPFMVXvs8HrcI%2FI8u%2FnAUMniX2a4RCxlZkSIYSSbVE2cVy8ETK08pHMLaE%2FciRr%2FUaGpeYJxDnCU%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241126T040419Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYZTSKULMN%2F20241126%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0d47411edac7df87f4b54d81da4665f3a317524af6db7eba11dcf6af4f9729a3&hash=6de1c4facec3d285a69f64febee99618e1bb28b9c69ec85e7261306f700654ee&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0196890409004695&tid=spdf-efeda5c7-5e71-44f6-970d-0408f2ae8751&sid=4afaa25e92dee841b98abc64158234a77a65gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f155d0a5b51090b505a00&rr=8e86fe65ead26a05&cc=us
  a0= 0.956; a1= 1.268; a2= 3.202; a3= -6.712; a4= 2.228; a5= -0.213; a6= 0.021; 
  # inner exponential expression
  inner_exp <- -exp(a2 + a3 * kt + a4 * kt^2 + a5 * m + a6 * m^2)
  # full function value
  k_value <- a0 - a1 * exp(inner_exp)
  return(k_value)
}

# simple airmass function (no dem)
airmass <- function(zenith){
  # m - airmass
  m = 1/cos(zenith * pi / 180)
  return(as.numeric(m) )
} 

# elevation-adjusted air mass based on Young (1987)
elevational_airmass <- function(zenith, elevation) {
  # base air mass calculation
  base_m <- 1 / cos(zenith * pi / 180)  # Convert zenith to radians
  # elevation adjustment using the barometric formula for pressure ratio
  pressure_ratio <- (1 - (2.25577e-5 * elevation))^5.2559
  # adjusted air mass
  adjusted_m <- base_m * pressure_ratio
  return(adjusted_m)
}

# calculate diffuse fraction (k0) using clearness, airmass, and diffuse fraction models
diffusef_Lk2_boulder <- function(IG, zenith_hr, I0 = 1367, m_elv = NULL) {
  # Calculate clearness index (kt)
  # kt <- clearness(IG, zenith_hr, I0)
  kt = IG/I0
  # Calculate airmass (m) 
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

# exo_irr(jdate)
k = diffusef_Lk2_boulder(IG = hrrr[['dswrf']], zenith, I0 = exo_irr(jdate), m_elv = hrrr[['gmp']])
plot(k)

plot(hrrr[['cld']])


k_modify_clouds <- function(k, Gh, cld_lyr){
  # Modify k for cloud layer
  Dh = 0 + 1*Gh
  
}

k_test = diffusef_Lk2_boulder(IG = sample(100:1000, 1000, replace = T), zenith, I0 = exo_irr(jdate))
kt_test = sample(seq(0.1,0.9, 0.01), 1000, replace = T)
m <- airmass(zenith)
k_test <- localG2(kt_test, m)
plot(k_test~kt_test)


kt = hrrr[['dswrf']]/exo_irr(jdate)


# RUN NEW GISPLIT4
k = diffusef_Lk2_boulder(IG = h[['dswrf']], zenith, I0 = exo_irr(jdate), m_elv = h[['gmp']])
newK = diffusef_gisplit4_cloudless(Gh = h[['dswrf']], 
                                   Ghc = sw3, Dhc = Id0, E0h = exo_irr(jdate)) 
plot(k)
plot(newK)

# Function to calculate K using the given equation
# After Ruiz-Arias and Gueymard, 2024 (GISPLIT)
diffusef_gisplit4_cloudless <- function(Gh, Ghc, Dhc, E0h) {
  # inputs
  # Gh; GHI or global horizontal irradiance 
  # Ghc: the modeled GHI under an ideal cloudless sky,  
  # Dhc: the modeled DIF under the same ideal cloudless sky 
  #
  # KT = Gh/E0h, where E0h is the extraterrestrial horizontal solar irradiance
  # Kcs = Gh/Ghc - clear sky index
  # Kds = Dhc/Ghc - clear sky diffuse fraction
  
  KT = Gh/E0h
  Kcs = Gh/Ghc
  Kds = Dhc/Ghc
  
  # Parameters for airmass -  f0; f1; f2; f3; f4
  fpars <- params_am_cloudless(m)
  
  # Calculate the argument inside the exponential
  exponent <- as.numeric(fpars['f0']) + as.numeric(fpars['f1']) * KT + as.numeric(fpars['f2']) * KT + as.numeric(fpars['f3']) * Kcs + as.numeric(fpars['f4']) * Kds
  
  # Calculate K using the logistic function
  K <- 1 / (1 + exp(exponent))
  
  return(K)
}

params_am_cloudless <- function(m){
  # tables in Ariaz-Ruiz and Gueymar 2024
  # https://pdf.sciencedirectassets.com/271459/1-s2.0-S0038092X24X0002X/1-s2.0-S0038092X24000574/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjENL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDlRw39ol1Pg5Ci%2FX5mIS2h%2BkV4%2BTTRezEEy%2BO8kqKYowIgIvAWA%2FJH8pW9o%2FnBwl5VSvrmBYlPt1LAYZZ6HGuMNugqswUIexAFGgwwNTkwMDM1NDY4NjUiDGAk8oYvCAdhtgemRSqQBZXh2ERfxbHH003wnxQx42e%2FfMW9Mp0vZ0%2FBitcOCIwkW20AxGUePDQxpmBj%2FAFMyvc6QfH1Iu5e6IXHUmZVFB4dotF0J5MvbMALV682am871zZ%2FTN18LQeLjfhcI%2FVujH9l5PmwM6eElff0q4UAWpr6cyadFJEI0idQSRGvpC8KewvOsg8H59mYbp44e7INXolX17y5jF83pM4IkBJ6ALbIrzwSaj%2B2QN5KPAfsjpR3Hs0OLf60Uqdwvu7zey4f35alJABgF9JEtSF%2FOEoJw3hb6VdYZz3HqkBVovgr9nfDDAIJ4QbNGhLqKwWcHB8tnKZ2m%2BBKkMJgxzNzn1apXuwjXpz8ExAd%2FFhDkN51k3GVNMzw85xiKdXxjyet3QF0VUYSm8WEg75Ft%2F9PCPxu%2FQ%2Fcy%2Bu4DlKZB2rc0L1hub13JBqvVe6H28RBUGKZeBGQKOZ87Czi0xgMl1BHnT4qcPP4yitnFOCnMmLm5w9VwzjV0NFYuRS%2BI3PGoyU%2BDXTghqgSCyyOACwCjQDMUJ%2FgA5Z%2BvV2%2FCUFwD01arcCAiSRCHgDGLRES9z3kcxLvQp7MfkMNc5eZnLhAhdQKdlSReHnOlKC7JPIVzEijwZYS5aIpsOatFzxo6CjcrHUE0LqAF5VOVFvhuVqJpjCqgyMygM3ZV5WsMRMu%2FfUfr%2FZ8BZlDZE86RQPq0xN%2BlSvbpo96wBQ%2B1msQKDnVSqALqiQSSOy6dHc9OAcqUX%2B8C55r%2F%2B9MqjhUDC7E8jo4T1l3LTv8U7EFMHpWkEgitmCycsQPPRjwkJb9k9s2bXeEfTnTH%2FCOgLAQHqDndlL6IGVMakiaRsO0M8wmg3t%2FWkeF6HwFx2WHayXhzdprSk%2FHHFOiP5mMMMaCqLoGOrEB9RRhQ2b2uXkEpRINHtzBohzafvIhgU63uZ9OkfqyV8NIwAZ2Ao5eFLCrwGqhxtki66wt45HfQflEH4VrOwhjwDtXjAd4%2Bq6ZnFD%2FMYy4n3m1Fgd5%2B4hvRan4d5yTWOsWg6zL1XpTXyZh75eUdz9zut87wnpsJaIVKaZiVD2IIo1l7su49lmVHtYbBbn192xkMULoLAO0G3VJn9mMhV0WGZaoZv3M%2FszzGzdLykdnj5yO&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241129T183348Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY7KZCLWXM%2F20241129%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e7a306d45318c8aecafec449f8f97edfde8f6c7eeb42a521d24b2828297a2f4d&hash=02eab031a706025fc9af0606f4051a8c904d434c78187eb7a7137b1c6d54ec27&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0038092X24000574&tid=spdf-e4c1d2bc-0729-47e2-b327-3d46784d55d7&sid=dc8d14b12521c24078197d9422b802d85440gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f155d0a02530d5e54095d&rr=8ea4b02f8ab327d4&cc=us
  am_vals1 = c(1.0,
               1.2115276586285884,
               1.4677992676220695,
               1.7782794100389228,
               2.154434690031884,
               2.610157215682537,
               3.1622776601683795,
               3.831186849557287,
               4.641588833612778,
               5.623413251903491,
               6.812920690579611,
               8.254041852680183)
  
  m_match = data.frame( range = 1:12, am_lower = am_vals1,
                        am_upper = c(am_vals1[2:12], 10.0) )
  
  # Create the list (equivalent to the provided Python list)
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
  
  # Convert the list to a matrix
  table_matrix <- do.call(rbind, table)
  
  # Convert the matrix to a data frame with column names 'am1' to 'am12'
  table_df <- as.data.frame(table_matrix)
  colnames(table_df) <- c("f0", "f1", "f2", "f3", "f4")
  table_df$am_range <- 1:12
  
  #
  m2 = m_match$range[m_match$am_lower <= m & m_match$am_upper > m]
  
  return(table_df[m2,1:5])
}



# Angle correction

# SASP 37°54’24.89088″N, -107°42’40.75924″W
cbind(37.906914, 107.711322) # https://snowstudies.org/swamp-angel-study-plot/
# SBB 37°54'25.30"N 107°42'40.23"W
# cbind(37.907028, -107.711175) # less accurate
cbind(37.906914, -107.726265) # https://snowstudies.org/senator-beck-study-plot/

angle2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  return(x)
}



# # python
# irwin_pyra = irwin_pyra.fillna({'VIS_IN': 0}).astype(np.float64)
# irwin_pyra[irwin_pyra < 0] = 0
# ## Max possible incoming from the sun
# irwin_pyra[irwin_pyra > 1300] = 1300
# ## no OUTCOMING when there is no INCOMING
# irwin_pyra.loc[irwin_pyra['VIS_IN'] <= 0, 'VIS_OUT'] = 0
# irwin_pyra.loc[irwin_pyra['NIR_IN'] <= 0, 'NIR_OUT'] = 0
# ## no OUTCOMING larger than INCOMING
# irwin_pyra['VIS_OUT'] = irwin_pyra.apply(lambda r: r['VIS_OUT'] if r['VIS_OUT'] < r['VIS_IN'] else r['VIS_IN'], axis=1)
# irwin_pyra['NIR_OUT'] = irwin_pyra.apply(lambda r: r['NIR_OUT'] if r['NIR_OUT'] < r['NIR_IN'] else r['NIR_IN'], axis=1)
# 
# irwin_pyra['net_solar'] = (
#   irwin_pyra['VIS_IN'] + irwin_pyra['NIR_IN'] - irwin_pyra['VIS_OUT'] - irwin_pyra['NIR_OUT']
# ).rolling(ROLL_WINDOW).mean().interpolate('time')
# 
# irwin_pyra['albedo'] = (
#   (0.67 * (irwin_pyra['VIS_OUT'] / irwin_pyra['VIS_IN'])) + (0.33 * (irwin_pyra['NIR_OUT'] / irwin_pyra['NIR_IN']))
# ).rolling(ROLL_WINDOW).mean().interpolate('time')
