#
# R instant insolation
#
# # # # # #

require(insol)
require(raster)
require(rgdal)


# # # # # # # # # # # # # # # # # # # #
#### TEST ############# TEST #########
# # # # # # # # # # # # # # # # # # # #

list.files("/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/")
list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/")

pth = list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes",
                 full.names = T)
file.exists(pth[522])

r = raster("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes/LC08_L2SP_035032_20150306_20200909_02_T1_SR_B3.TIF")
plot(r)

rgb1 = list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes/", 
                  pattern="035032_20150306_.*_B[1-9].TIF", full.names = T)
rgb2 = stack(rgb1)
plotRGB(rgb2, 3,2,1, stretch = "lin")

# create buffer point - viz (TEST)
# pts = click(n=1)
pts = cbind(301796.9, 4437015)
pts2 <- SpatialPoints(coords = pts, proj4string = CRS(proj4string(rgb2)) )
plot(pts2, add=T,col='red', cex=3, pch=20)
pb = buffer(pts2, 1.6e4)
plot(extent(pb), add=T, col='red')

rc = crop(rgb2, pb)
plotRGB(rc, 3,2,1, stretch = "lin")

# DEM - read (CROP)
d =  raster("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/fusion_files/DEM/034033_merged.tif")
plot(d)

pb_proj = spTransform(pb, proj4string(d) )
d2 = crop(d, pb_proj)
d2 = projectRaster(d2, crs= proj4string(rgb2))
plot(d2)

# check transformation in subregion
dpoly = as(extent(d2), 'SpatialPolygons')
crs(dpoly) = crs(d2)
plot(r)
plot(dpoly, add=T)



# # # # # # # # # # # # # # # # # # # #
#### TEST ############# TEST #########
# # # # # # # # # # # # # # # # # # # #


### WORKFLOW (for single subset) #####
# 1. For each tile 
# 2. extract all vars
# 3. compute moment insol

list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes",
                 full.names = T, pattern = "MTL.txt")
## METADATA
mtl_path = "/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes/LC08_L2SP_035032_20150306_20200909_02_T1_MTL.txt"
con = file(mtl_path)
line = readLines(con)
#
## SOLAR POSITION AND DATETIME
ll = trimws(grep("SUN_ELEVATION", line, value = T))
s_elv = as.numeric(strsplit(ll, "= ")[[1]][2])
s_elv
s_az = as.numeric(strsplit(trimws(grep("SUN_AZ", line, value = T)), "= ")[[1]][2])
s_az
# ACQUISITION TIME
# trimws(grep("DATE_ACQUI", line, value = T))
# trimws(grep("TIME", line, value = T))
date_time = paste(strsplit(trimws(grep("DATE_ACQUI", line, value = T)), "= ")[[1]][2],
                   substr(strsplit(strsplit(trimws(grep("TIME", line, value = T)), "= ")[[1]][2], "\"")[[1]][2], 1, 8))
dtime = as.POSIXct(date_time, tz = "UTC")
#
## REQUIRED TERRAIN LAYERS
p0 <<- as(extent( c(294458.7, 326458.7, 4411268, 4443268) ), 'SpatialPolygons') # same as d2
crs(p0) <- "+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs"
if(is.null(d)){stop()} # check dem loaded 
if(is.null(d2)){stop()} else{
  d2 = crop(d, pb_proj)
  d2 = projectRaster(d2, crs= proj4string(rgb2))
}
d2
d_mat <- as.matrix(d2)
lat_lon <- c(round((d@extent@ymax + d@extent@ymin)/2,5), round((d@extent@xmax + d@extent@xmin)/2,5))
dem_res <- dem.res(d2, lat_lon[1]) # if in lat
dem_res <- mean(res(d2))
VF_mat <- view.factor(d_mat, d2, dem_res)
tmz <- round(lat_lon[2]/15,1)
s_a <- slope.aspect(d2)

# shade, insol, illumination

## SHADE
sol_vect <- insol::normalvector(90-s_elv,s_az)
# sol_vect2 <-  c(-sol_vect[1:2],sol_vect[3]))

shd <- insol::doshade(d2,sol_vect)
# sh_ls <- flip(t(flip(t(sh_ls),1)),1)

# test
plot(d2)
plot(shd, add=T, col=c("green",NA))
im_crop = crop(rgb2, d2)
plotRGB(im_crop, 3,2,1)
plot(shd, add=T, col=c("green",NA))

# edit shade
# shd2 = shd*!is.na(ld)
shd2 = shd;shd2[is.na(ld)] = NA
# shd2 = !shd2; shd2[shd2==0] = NA
plot(shd2)
#
# filter pixels based on 'shade_size_filter'
# size etc
shade_size_filter = 20; window_size=9
shclump <-  clump(!shd2, directions=8)
f<-as.data.frame(freq(shclump))
exludeShade <- f$value[which(f$count <= shade_size_filter)]
shfilter <- shclump
shfilter[shclump %in% exludeShade] <- NA 
shf <- (!is.na(shfilter))*!is.na(d2)
plot(shf)
# xompare
plot(shd)
plot(shf, add=T, col=c(NA,'red'))
plotRGB(im_crop, 3,2,1)
plot(shd, add=T, col=c("green",NA))
plot(shf, add=T, col=c(NA,'red'))

# final shade
!shf


#### TEST full day
# zenith and azimuth angles
yr = format(dtime, "%Y")
mo = format(dtime, "%m")
dy = format(dtime, "%d")
jd=JD(seq(ISOdate(yr,mo,dy,0, tz = "MST"),ISOdate(yr,mo,dy,23, tz="MST"),by="1 hour"))
jd=JD(seq(ISOdate(yr,mo,dy,0),ISOdate(yr,mo,dy,23),by="1 hour"))

# jd = JD(dtime)
JD(jd, inverse=T)

# sun position and vector
sv = sunvector(jd,lat_lon[1],lat_lon[2],-7); sp1=sunpos(sv)
# daylight hours (zenith <= 90)
sp=sp1[which(sp1[,2]<=90),]
sv=sv[which(sp1[,2]<=90),]

zenith=sp[,2]
az_noon = which.min(abs(180-sunpos(sv))[,1])
azimuth_eq = c(sunpos(sv)[,1][1:az_noon]-180,sunpos(sv)[,1][(az_noon+1):length(sunpos(sv)[,1])]-180)

### for moment
i = 5

# zenith and incident angles
cos_inc <- cos.slope(zenith[i], azimuth_eq[i], aspect = s_a[[2]], slope = s_a[[1]])
cos_sfc <- cos(radians(zenith[i]))
plot(make.raster(cos_inc, d2) )


# test
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100) )
for (i in 1:length(zenith)){
  cos_inc <- cos.slope(zenith[i], azimuth_eq[i], aspect = s_a[[2]], slope = s_a[[1]])
  sh0 <- doshade(d_mat, sv[i,], dl=30)
  plot(hsh, col=grey(c(0:100)/100), main = paste("zenith:", round(zenith[i]), "| time:", JD(jd[i], inverse=T)), legend=F)
  plot(make.raster(cos_inc, d2), add = T, col = terrain.colors(100, 0.4))
  plot(make.raster(sh0, d2), add=T, col=c("black",NA), legend=F)
}

# single moment:
## VARIABLES
zenith0 = 90-s_elv
azimuth_eq0 = s_az - 180
cos_inc0 <- cos.slope(zenith0, azimuth_eq0, aspect = s_a[[2]], slope = s_a[[1]])
jd0 = JD(dtime)
height <- array(d2)
# height <- cellStats(d2, 'mean')
visibility <- 28 
RH <- 60
tempK <- 278.15 # 5 C or 41 F

# ! ! ! ! ! ! ! ! ! ! ! 
# INSOLATION #### !!!!!
# insolation arriving perpendicular to solar beam (direct and diffuse)
Idirdif = insolation(zenith0,jd0,height,visibility,RH,tempK,0.002,0.35)
Ib = matrix(Idirdif[,1],nrow=nrow(d2),ncol=ncol(d2), byrow = T)
Id = matrix(Idirdif[,2],nrow=nrow(d2),ncol=ncol(d2), byrow=T)

# Ib0 = make.raster(Ib * cos_inc0, d2)
cosi = make.raster(cos_inc0,d2)
Ib0 = make.raster(Ib, d2)
Id0 = make.raster(Id, d2) 

# swin = (Ib * shf + Id)*cos_inc
# (Ib * sh + Id * VF_mat )*cos_inc
sw1 = Ib0 * cosi * !shf # direct beam only
sw2 = sw1 + Id0 # direct and diffuse
sw3 = sw1 + (Id0 * make.raster(VF_mat, d2)) # includes view factor

# final insol plot
# plot(make.raster(VF_mat, d2), col=grey(c(0:100)/100) )
# plot(hsh, col=grey(c(0:100)/100), main = paste("zenith:", round(zenith[i]), "| time:", JD(jd[i], inverse=T)), legend=F)
# plot(make.raster(cos_inc, d2), add = T, col = terrain.colors(100, 0.4))

# model 1
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw1, col=colorspace::heat_hcl(100, alpha=0.5), add=T, axes=FALSE, 
     legend.args = list(text = 'Insolation (Wm-2)', side = 3, 
                        font = 2,  line = 0, cex = 0.8), mar = c(2, 3, 3, 3), )

# model 2
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw2, add = T, col = heat.colors(100, 0.5))

# model 3
par()$mar
par(mar = c(1,1,1,4))
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw3, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
r.range <- c(minValue(sw3), maxValue(sw3))
plot(sw3, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=seq(r.range[1], r.range[2], 100),
                    labels=round(seq(r.range[1], r.range[2], length.out = 10)), 
                    cex.axis=0.6),
     legend.args=list(text='Insolation (Wm^2)', side=3, font=2, line=2.5, cex=1))



# FINAL PLOT
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw1, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
r.range <- c(minValue(sw3), maxValue(sw3))
plot(sw1, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
legend("bottomleft", c("Clear-sky Modeled Insolation",
                       paste("TilePath:", strsplit(basename(mtl_path), "_")[[1]][3], "| Date:",
                             format(dtime, tz="MST",usetz=TRUE) ) ) )
# FINAL PLOT
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw3, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
r.range <- c(minValue(sw3), maxValue(sw3))
plot(sw3, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
legend("bottomleft", c("Clear-sky Modeled Insolation",
                       paste("TilePath:", strsplit(basename(mtl_path), "_")[[1]][3], "| Date:",
                             format(dtime, tz="MST",usetz=TRUE) ) ) )


# time
format(dtime, "%H:%M:%S")
format(dtime, tz="MST",usetz=TRUE)


### FUNCTIONS



# make extract vars function
extract.solar.vars <- function(ls_path, d_path, tilepath, subset0 = NULL){
  #
  # complies DEM, LS RGB
  #
  # get variables
  
  # lns = read.delim(file_path) # gsub("([0-9]+)\\..*$", "\\1", ll)
  con = file(mtl_path)
  line = readLines(con)
  # SUN_ELEVATION
  ll = trimws(grep("SUN_ELEVATION", line, value = T))
  s_elv = as.numeric(strsplit(ll, "= ")[[1]][2])
  s_elv
  # ACQUISITION TIME
  trimws(grep("DATE_ACQUI", line, value = T))
  trimws(grep("TIME", line, value = T))
  trimws(grep("SUN_AZ", line, value = T))
  # trimws(grep("ZEN", line, value = T))
  
  # make extract sun angles function
  
}

###############################################
slope.aspect <- function(dem, units = 'radians', neighbor = 8){
  s <- terrain(dem, opt='slope',unit=units,neighbors=neighbor)
  a <- terrain(dem, opt='aspect',unit=units,neighbors=neighbor)
  stk <- stack(s,a)
  return(stk)
}

###############################################
cos.slope <- function(zenith, azimuth_eq, aspect, slope){
  # returns a matrix of the cosine of the incident angle at a given moment
  exposures  = aspect - radians(180)
  cos_inc = acos((cos(radians(zenith)) * cos(slope)) +
                   (sin(radians(zenith)) * sin(slope) * cos(radians(azimuth_eq) - exposures)))
  
  cos_inc = as.matrix(cos_inc)
  # get rid of self shading values
  cos_inc[cos_inc > radians(90)] = radians(90)
  cos_inc = cos(cos_inc)
  return(cos_inc)
}

###############################################
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}

###############################################
view.factor <- function(dem_mat, dem, dem_res, elv_interval = 5, az_interval = 15){
  print("____Generating sky view factor_____")
  ptm <- proc.time()
  ELV = rev(seq(0, 90, elv_interval))
  AZI = seq(0, 345, az_interval)
  AZ = matrix(0,nrow=dim(dem_mat)[1]*dim(dem_mat)[2],ncol=length(AZI))
  for (vv in 1:length(AZI)){
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

###############################################
dem.res <- function(dem, latitude){
  resolution <- res(dem)[1]/(1/(111320*cos(radians(latitude))))
  dem_res <- round(resolution, 2)
  return(dem_res)
}

# make shade function

# make insol instant function
