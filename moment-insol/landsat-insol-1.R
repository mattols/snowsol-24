#
# R instant isolation
# FUNCTIONS
# 
# # # # # #

require(insol)
require(raster)
require(rgdal)
require(rgeos)

# # # # # PATHS
main_path = "/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes"
# save_path = "/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/drfs/oli-drfs/agu24/insol"
save_path = "/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/drfs/oli-drfs/agu24/insol2"
dpath_main = "/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/fusion_files/DEM"



###############################################
### FUNCTIONS

###############################################
sw.solar.noon <- function(mtl_path, save_sw = TRUE){
  # SOLAR NOON
  # zenith and azimuth angles
  yr = format(dtime, "%Y"); mo = format(dtime, "%m"); dy = format(dtime, "%d")
  jd_mst=JD(seq(ISOdate(yr,mo,dy,0, tz = "MST"),ISOdate(yr,mo,dy,23, tz="MST"),by="1 mins"))
  jd=JD(seq(ISOdate(yr,mo,dy,0),ISOdate(yr,mo,dy,23),by="1 mins"))
  # sun position and vector
  # timezone
  zchng = dtime; attr(zchng, "tzone") <- "US/Mountain"
  tmz = as.numeric(format(zchng, "%H")) - as.numeric(format(dtime, "%H"))
  #
  sv = sunvector(jd,lat_lon[1],lat_lon[2],tmz); sp1=sunpos(sv)
  # daylight hours (zenith <= 90)
  sp=sp1[which(sp1[,2]<=90),]
  sv=sv[which(sp1[,2]<=90),]
  # solar noon
  az_noon = which.min(abs(180-sunpos(sv))[,1])
  zenith = sp[az_noon,2]
  azim = sp[az_noon,1]
  azimuth_eq = azim -180
  
  jd0 = jd[which(sp1[,2]<=90)][az_noon]
  
  # # provide moment?
  # lz_over = which.min(abs((90-s_elv)-sunpos(sv))[,2])
  # la_over = which.min(abs((s_az)-sunpos(sv))[,1])
  # zenith = sp[lz_over,2]
  # azim = sp[la_over,1]
  
  # rerun shade
  # make.shade(zenith, azim)
  if(dim(d)[1] != dim(d_mat)[1]){d_mat <<- as.matrix(d)}
  sh0 <- doshade(d_mat, sv[az_noon,], dl=30)
  sh0 <- make.raster(sh0, d) 
  # zenith and incident angles
  cos_inc <- cos.slope(zenith, azimuth_eq, aspect = s_a[[2]], slope = s_a[[1]])
  
  # SOLAR NOON INSOLATION 
  height <- array(d)
  Idirdif = insolation(zenith,jd0,height,visibility,RH,tempK,O3,alphag)
  Ib = matrix(Idirdif[,1],nrow=nrow(d),ncol=ncol(d), byrow = T)
  Id = matrix(Idirdif[,2],nrow=nrow(d),ncol=ncol(d), byrow=T)
  # MAKE RASTER
  # cosi = make.raster(cos_inc0,d)
  Ib0 = make.raster(Ib, d)
  Id0 = make.raster(Id, d) 
  
  sw1 = Ib0 * cos_inc * sh0 # direct beam only
  sw2 = sw1 + Id0 * cos(radians(zenith)) # direct and diffuse
  # sw3 = sw1 + (Id0 * make.raster(VF_mat, d)) # includes view factor
  
  # # proj
  # pat <- substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_path))), 1, 15)
  # rpath <- list.files(main_path, full.names = T, pattern = paste0(pat,".*B[1-7].TIF"))[1]
  # r <- raster(rpath)
  # sw2 <-  projectRaster(sw2, crs= proj4string(r))
  
  # time 
  dmoment = JD(jd_mst[which(sp1[,2]<=90)], inverse=T)[az_noon]
  dmoment_char_noon <<- gsub(":", "-", gsub("\\s","_", as.character(dmoment) ))
  
  if(isTRUE(save_sw)){
    sname = substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_path))), 1, 15)
    writeRaster(sw2, file.path(save_path, paste0("dsw_", sname, "_solnoon_", dmoment_char_noon, ".tif")) )
  }else{return(sw2)}
  
}

###############################################
insol.params <- function(){
  visibility <<- 25 
  RH <<- 40 # 60
  tempK <<- 273.15 #278.15 # 5 C or 41 F
  O3 <<- 340/1000 # DU = 0.01 mm (CO 325-340 DU)
  alphag <<- 0.35
  # https://disc.gsfc.nasa.gov/datasets?page=1&source=AURA%20OMI
}

###############################################
sw.moment.mtl <- function(mtl_path, save_sw = TRUE){
  ## VARIABLES
  zenith0 = 90-s_elv
  azimuth_eq0 = s_az - 180
  cos_inc0 <- cos.slope(zenith0, azimuth_eq0, aspect = s_a[[2]], slope = s_a[[1]])
  jd0 = JD(dtime)
  # ! ! ! ! ! ! ! ! ! ! ! 
  # INSOLATION
  # insolation arriving perpendicular to solar beam (direct and diffuse)
  height <- array(d)
  Idirdif = insolation(zenith0,jd0,height,visibility,RH,tempK,O3,alphag)
  Ib = matrix(Idirdif[,1],nrow=nrow(d),ncol=ncol(d), byrow = T)
  Id = matrix(Idirdif[,2],nrow=nrow(d),ncol=ncol(d), byrow=T)
  # MAKE RASTER
  # cosi = make.raster(cos_inc0,d)
  Ib0 = make.raster(Ib, d)
  Id0 = make.raster(Id, d) 
  
  # make shade
  if(!exists("shdd")){make.shade2(90-s_elv, s_az)}
  
  # swin = (Ib * shf + Id)*cos_inc
  # (Ib * sh + Id * VF_mat )*cos_inc
  sw1 = Ib0 * cos_inc0 * shdd # direct beam only
  sw2 = sw1 + Id0 * cos(radians(zenith0)) # direct and diffuse
  # sw3 = sw1 + (Id0 * make.raster(VF_mat, d)) # includes view factor
  
  dmoment_char <<- gsub(":", "-", gsub("\\s","_", as.character(dtime) ))
  
  if(isTRUE(save_sw)){
    sname = substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_path))), 1, 15)
    writeRaster(sw2, file.path(save_path, paste0("dsw_", sname, "_moment_", dmoment_char, ".tif")) )
  }else{return(sw2)}
}

###############################################
make.shade <- function(ZEN, AZ, shade_size_filter = 40){
  ## SHADE
  sol_vect <- insol::normalvector(ZEN,AZ)
  shd <- insol::doshade(d,sol_vect)
  shd2 = shd;shd2[is.na(d)] = NA
  #
  # filter pixels based on size
  shclump <-  clump(!shd2, directions=8)
  f<-as.data.frame(freq(shclump))
  exludeShade <- f$value[which(f$count <= shade_size_filter)]
  shfilter <- shclump
  shfilter[shclump %in% exludeShade] <- NA 
  shf <- (!is.na(shfilter))*!is.na(d)
  shdd <<- !shf
}

###############################################
make.shade2 <- function(ZEN, AZ){
  ## SHADE
  # ZEN = 90 - s_elv ; AZ = s_az
  sol_vect <- insol::normalvector(ZEN,AZ)
  # find matching raster 
  pat <- substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_path))), 1, 15)
  rpath <- list.files(main_path, full.names = T, pattern = paste0(pat,".*B[1-7].TIF"))[1]
  r <- raster(rpath)
  d <<-  projectRaster(d, crs= proj4string(r))
  s_a <<-  projectRaster(s_a, crs= proj4string(r))
  shdd <<- insol::doshade(d, sol_vect)
}


###############################################
extract.solar.vars <- function(mtl_path){
  # GET SOLAR VARS
  ## METADATA
  con = file(mtl_path)
  line = readLines(con)
  #
  ## SOLAR POSITION AND DATETIME
  ll = trimws(grep("SUN_ELEVATION", line, value = T))
  s_elv <<- as.numeric(strsplit(ll, "= ")[[1]][2])
  s_az <<- as.numeric(strsplit(trimws(grep("SUN_AZ", line, value = T)), "= ")[[1]][2])
  # ACQUISITION TIME
  # trimws(grep("DATE_ACQUI", line, value = T))
  # trimws(grep("TIME", line, value = T))
  date_time <<- paste(strsplit(trimws(grep("DATE_ACQUI", line, value = T)), "= ")[[1]][2],
                    substr(strsplit(strsplit(trimws(grep("TIME", line, value = T)), "= ")[[1]][2], "\"")[[1]][2], 1, 8))
  dtime <<- as.POSIXct(date_time, tz = "UTC")
}

###############################################
terrain.variables <- function(mtl_path, dpath_main, minimize=TRUE){
  #
  # create required terrain data for insol
  # minimize crops to smallest possible area of overlap
  #
  pathrow00 <- strsplit(basename(mtl_path), "_")[[1]][3]
  # if (pathrow00=="034032") {pathrow00 = "034033"}
  # PATHS
  # ppprrr <- strsplit(basename(mtl_path), "_")[[1]][3]
  # dfiles <- list.files(dpath_main, full.names = T, pattern = "merged.tif")
  # dpath <- dfiles[pmatch(pathrow0, basename(dfiles))]
  #
  dfiles <- list.files(dpath_main, full.names = T, pattern = pathrow00, recursive = T)
  # asp_dem_slop <-  dfiles[grepl("merged|aspect|slope", basename(dfiles))]
  asp_dem_slop <-  dfiles[grepl("merged", basename(dfiles))]
  # # save_name <- file.path(save_path, paste(pathrow0, "_terrain.tif"))
  # if (file.exists(save_name)){
  #   print("terrain layers exist")
  #   # READ FILE
  # }else{
  #   print("Generating terrain layers...")
  #   # OPEN
  #   # TERRAIN
  #   # SAVE FILE
  # }
  d <<- raster(asp_dem_slop[[1]])
  # d2 <- projectRaster(d2, crs= proj4string(rgb2))
  lat_lon <<- c(round((d@extent@ymax + d@extent@ymin)/2,5), round((d@extent@xmax + d@extent@xmin)/2,5))
  # resolution
  if (trimws(strsplit(as.character(crs(d)), "\\+")[[1]][2]) == "proj=longlat"){
    dem_res <- dem.res(d, lat_lon[1]) # if in lat
  }else{dem_res <- mean(res(d)) } # if UTM
  dem_res <<- 30 # manually set
  
  # VF_mat <- view.factor(d_mat, d, dem_res)
  # tmz <- round(lat_lon[2]/15,1)
  # # timezone
  # zchng = dtime; attr(zchng, "tzone") <- "US/Mountain"
  # tmz = as.numeric(format(zchng, "%H")) - as.numeric(format(dtime, "%H"))
  # s_a <<- stack(asp_dem_slop[[3]], asp_dem_slop[[1]])
  # writeRaster(stack(stack(d, s_a), make.raster(VF_mat, d), file.path(save_path, paste(pathrow0, "_terrain.tif"))))
  if(minimize){
    if(!exists("opoly")){
      overlap.ls.dem(mtl_path, dmain_path)
      opoly2 <- spTransform(opoly, proj4string(d) )
      d <<- crop(d, opoly2)
      # s_a <<- crop(s_a, opoly2)
    }
  }
  s_a <<- slope.aspect(d)
  d_mat <<- as.matrix(d)
  lat_lon <<- c(round((d@extent@ymax + d@extent@ymin)/2,5), round((d@extent@xmax + d@extent@xmin)/2,5))
  
}

###############################################
  ###############################################
############################################### complete
  ###############################################


###############################################
list.scene.info <- function(main_path = "/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes"){
  mtl_files <<-  list.files(main_path, full.names = T, pattern = "MTL.txt")
  # scene info
  pathrows <<- unique(unlist(lapply(1:length(mtl_files), function(x) strsplit(basename(mtl_files[x]), "_")[[1]][3])))
  dates0 <<- unique(unlist(lapply(1:length(mtl_files), function(x) strsplit(basename(mtl_files[x]), "_")[[1]][4])))
  all_dates <<- substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_files))), 1, 15)
  print(paste(length(pathrows), "unique scene locations |", length(all_dates), "total images"))
  
  # limit based on existing DEMs
  # return(mtl_files)
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
cos.slope2 <- function(zenith, azimuth_eq, aspect, slope, as_mat=FALSE){
  # simplified does not work as well
  exposures  = aspect - radians(180)
  cos_inc = sin(radians(zenith)) * sin(slope) +
                   cos(radians(zenith)) * cos(slope) * cos(radians(azimuth_eq) - exposures)
  
  if(as_mat){cos_inc = as.matrix(cos_inc)}
  # get rid of self shading values
  cos_inc[acos(cos_inc) > radians(90)] = radians(90)
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

###############################################
dem.res <- function(dem, latitude){
  resolution <- res(dem)[1]/(1/(111320*cos(radians(latitude))))
  dem_res <- round(resolution, 2)
  return(dem_res)
}

# minimize computation
###############################################
overlap.ls.dem <- function(mtl_path, dmain_path){
  ##
  # determine minimum extent for each scene - extract and save overlap
  #
  # find matching raster 
  pat <- substr(unique(gsub(".*L2SP_(.*)_02.*", "\\1", basename(mtl_path))), 1, 15)
  rpath <- list.files(main_path, full.names = T, pattern = paste0(pat,".*B[1-7].TIF"))[1]
  r <- raster(rpath)
  rpoly = as(extent(r), 'SpatialPolygons')
  crs(rpoly) = crs(r)
  # d2 <-  projectRaster(d, crs= proj4string(r))
  # gather all dimentions for DEMs
  if(!exists("d")){
    d = raster(list.files(dpath_main, full.names = T)[1])
  }
  # polygonize
  dpoly = as(extent(d), 'SpatialPolygons')
  crs(dpoly) = crs(d)
  dpoly_proj = spTransform(dpoly, proj4string(r) )
  
  # overlap area
  opoly <<- rgeos::gIntersection(dpoly_proj,rpoly)
  
  # test
  # plot(rpoly)
  # plot(dpoly_proj,border='red',add=T)
  # plot(rgeos::gIntersection(dpoly_proj,rpoly), add=T,col='blue',lwd=2)
}









# NOT NEEDED
overlap.roi.pathrow <- function(pathrows, mtl_path, dmain_path){
  ##
  # determine for each scene - extract and save overlap
  #
  
  # gather all dimentions for DEMs
  d = raster(list.files(dpath_main, full.names = T)[1])
  dpoly = as(extent(d), 'SpatialPolygons')
  crs(dpoly) = crs(d)
  dpoly_proj = spTransform(dpoly, proj4string(r) )
  
  for (i in 1:length(pathrows)){
    i = 2
    # take first
    r = raster(list.files(main_path, full.names = T, pattern = paste0(pathrows[i],".*B[1-7].TIF"))[1])
    # r2 = !is.na(r);r2[r2==0] = NA
    # rpoly = rasterToPolygons(r2, dissolve=TRUE)
    rpoly = as(extent(r), 'SpatialPolygons')
    crs(rpoly) = crs(r)
    opoly <- rgeos::gIntersection(dpoly_proj,rpoly)
    
    pathrows <<- unique(unlist(lapply(1:length(mtl_files), function(x) strsplit(basename(mtl_files[x]), "_")[[1]][3])))
  }
  
  plot(rpoly)
  plot(dpoly_proj,border='red',add=T)
  plot(rgeos::gIntersection(dpoly_proj,rpoly), add=T,col='blue',lwd=2)
  
}

