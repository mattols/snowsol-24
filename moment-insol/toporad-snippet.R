# TOPORAD - functions

# 
# TOPOGRAPHIC SOLAR RADIATION MODEL FUNCTIONS
# 01/28/2019
# 
# Matthew Olson - University of Utah, Department of Geography 
# email: matthew.olson@geog.utah.edu
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

require(insol)
require(raster)
require(rgdal)

list.files("/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/isnobal-data/")

pth = list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes",
                 full.names = T)
file.exists(pth[522])

r = raster("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes/LC08_L2SP_035032_20150306_20200909_02_T1_SR_B3.TIF")
plot(r)

rgb1 = list.files("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/OLI_snow_scenes/", 
                  pattern="035032_20150306_.*_B[1-9].TIF", full.names = T)
rgb2 = stack(rgb1)
plotRGB(rgb2, 3,2,1, stretch = "lin")

# create buffer point - viz
pts = click(n=1)
# pts = cbind(301796.9, 4437015)
pts2 <- SpatialPoints(coords = pts, proj4string = CRS(proj4string(rgb2)) )
plot(pts2, add=T,col='red', cex=3, pch=20)
pb = buffer(pts2, 1.6e4)
plot(extent(pb), add=T, col='red')

rc = crop(rgb2, pb)
plotRGB(rc, 3,2,1, stretch = "lin")

d =  raster("/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DRFS-Fusion/fusion_files/DEM/034033_merged.tif")
plot(d)

pb_proj = spTransform(pb, proj4string(d) )
d2 = crop(d, pb_proj)
d2 = projectRaster(d2, crs= proj4string(rgb2))
plot(d2)

# required layers
# DEM
demtype="NASADEM"
demr <- raster(paste0("/uufs/chpc.utah.edu/common/home/cryosphere/molson/DebrisCover/data/Landsat/resultSR/140041/140041","_dem",demtype,"30.grd"))
p2 <<- as(extent( c(457000, 490000, 3088000, 3120000) ), 'SpatialPolygons')
crs(p2) <- "+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs";p2<<-p2
dem = crop(demr, p2)
#
# make globally available
d_mat <<- as.matrix(dem)
lat_lon <<- c(round((dem@extent@ymax + dem@extent@ymin)/2,5), round((dem@extent@xmax + dem@extent@xmin)/2,5))
dem_res <<- dem.res(dem, lat_lon[1])
VF_mat <<- view.factor(d_mat, dem, dem_res)
tmz <<- round(lat_lon[2]/15,1)
s_a <<- slope.aspect(dem)


topo.shade <- function(tile_list, save_name,demtype=c("NASADEM","SRTM")){
  save_base_name <- strsplit(save_name,"_")[[1]][1]
  ble  <- tile.raster.list(tile_list,"B2")
  # minimizing to glaciers
  gp2 <- fasterize::fasterize(readRDS(paste0(save_base_name,"_glacierPolygons.rds")),ble[[1]]) %>% mask(ble[[1]]) %>% trim()
  ble2 <- lapply(ble, function(x) x %>% crop(.,gp2) )
  # get dem if doesn't exist
  if(!file.exists(paste0(save_base_name,"_dem",demtype,"30.grd"))){
    if(!file.exists(paste0(save_base_name,"_glacierPolygons.rds"))){stop("Must run `g.poly()` to determine DEM extent")}
    find.dem(save_base_name, type=demtype)
  }
  demr <- raster(paste0(save_base_name,"_dem",demtype,"30.grd")) %>% crop(ble2[[1]])
  dem_res <- lapply(ble2, function(y) demr%>%resample(y))
  d_flip <- lapply(dem_res, function(d) flip(t(raster::flip(t(d),1)),1))
  # solar geom
  xm_ls <- lapply(grep("MTL.xml",list.files(tile_list, full.names = T),value=TRUE), 
                  function(m) m %>% xml2::read_xml() )
  selv_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_ELEVATION")%>%xml2::xml_text()%>%as.numeric() )
  saz_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_AZIMUTH")%>%xml2::xml_text()%>%as.numeric() )
  sv_ls <- lapply(1:length(selv_ls),function(i) insol::normalvector(90-selv_ls[[i]][1],saz_ls[[i]][1]))
  sv_ls_cor <- lapply(sv_ls, function(sv) c(-sv[1:2],sv[3]))
  sh_ls <- lapply(1:length(d_flip), function(s) insol::doshade(d_flip[[s]],sv_ls_cor[[s]]))
  sh_ls <- lapply(sh_ls, function(s) flip(t(flip(t(s),1)),1))
  # extend to full tile
  sh_ls <- lapply(sh_ls, function(x) extend(x, extent(ble[[1]])))
  return(do.call(brick,sh_ls))
}

###############################################
dem.res <- function(dem, latitude){
  resolution <- res(dem)[1]/(1/(111320*cos(radians(latitude))))
  dem_res <- round(resolution, 2)
  return(dem_res)
}

###############################################
slope.aspect <- function(dem, units = 'radians', neighbor = 8){
  s <- terrain(dem, opt='slope',unit=units,neighbors=neighbor)
  a <- terrain(dem, opt='aspect',unit=units,neighbors=neighbor)
  stk <- stack(s,a)
  return(stk)
}

###############################################
location.variables <- function(demL, shape, resampleFactor = 1){
  # varibles that remain constant over time
  # crop raster
  dem <<- crop.raster(demL, shape)
  dem <<- void.fill(dem)
  glacier <<- shape
  # resample if necessary
  if (resampleFactor == 1){
    print("Keep native resolution")
  } else{
    dem <<- new.resolution(dem, resampleFactor)
  }
  # make globally available
  d_mat <<- as.matrix(dem)
  lat_lon <<- c(round((dem@extent@ymax + dem@extent@ymin)/2,5), round((dem@extent@xmax + dem@extent@xmin)/2,5))
  dem_res <<- dem.res(dem, lat_lon[1])
  VF_mat <<- view.factor(d_mat, dem, dem_res)
  tmz <<- round(lat_lon[2]/15,1)
  s_a <<- slope.aspect(dem)
  print("Location variables loaded to memory")
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
create.tfmodels <- function(){
  # generate empty model matrices
  # model_flat      <<- array(0,dim=dim(d_mat))
  # model_flat_base <<- array(0,dim=dim(d_mat))
  # model_flat_sh   <<- array(0,dim=dim(d_mat))
  # model_inc_sr    <<- array(0,dim=dim(d_mat))
  # model_inc_sh    <<- array(0,dim=dim(d_mat))
  model_vf_base   <<- array(0,dim=dim(d_mat))
  model_vf        <<- array(0,dim=dim(d_mat))
  model_ref       <<- array(0,dim=dim(d_mat)) 
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
atmos.vars <- function(){
  # atmospheric variables (constant simulating clear-sky)
  height <<- array(d_mat)
  visibility <<- 28 
  RH <<- 60
  tempK <<- 278.15
}

###############################################
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}

###############################################
crop.raster <- function(stk, shp, buffer = 0.03){
  r <- crop(stk,extent(xmin(shp)-buffer,xmax(shp)+buffer,
                       ymin(shp)-buffer,ymax(shp)+buffer))
  return(r)
}




# # # # # # # # # # 
 # # # # # # # # # # 
# # # # # # # # # # 

   
swtest = sw.daily(sw_totals = T)

###############################################
sw.daily <- function(date = ISOdate(2017, 6, 21, 0), perc = FALSE, sw_totals = FALSE, plot_moment = FALSE){
  # daily sw variables -> call sw.moment()
  year  <- format(date,'%Y')
  month <- format(date,'%m')
  day   <- format(date,'%d')
  print(paste("Calculating shortwave for", ISOdate(year, month, day, 0)))
  ptm <- proc.time()
  
  # julian moment at every 15 mins
  # jd=JD(seq(ISOdate(year,month,day,0),ISOdate(year,month,day,23),by="15 mins"))
  jd=JD(seq(ISOdate(year,month,day,0),ISOdate(year,month,day,23),by="1 hour"))
  
  # sun position and vector
  sv = sunvector(jd,lat_lon[1],lat_lon[2],5.5); sp1=sunpos(sv)
  
  # daylight hours (zenith <= 90)
  sp=sp1[which(sp1[,2]<=90),]
  sv=sv[which(sp1[,2]<=90),]
  
  # zenith and azimuth angles
  zenith=sp[,2]
  az_noon = which.min(abs(180-sunpos(sv))[,1])
  azimuth_eq = c(sunpos(sv)[,1][1:az_noon]-180,sunpos(sv)[,1][(az_noon+1):length(sunpos(sv)[,1])]-180)
  
  ## LOOP
  atmos.vars()
  create.tfmodels()
  print(paste("...cycling through", length(zenith), "moments"))
  for (m in 1:length(zenith)){
    # calculate fluxes at time m
    create.tfmodels()
    sw.moment(sv[m,], zenith[m], azimuth_eq[m], jd[m])
    if (m==1){stkmod = make.raster(model_ref, dem)}else{
      stkmod = stack(stkmod, make.raster(model_ref, dem))
    }
    plot(stkmod[[m]], col=heat.colors(40))
    # sh0 = make.raster(sh, d_flip)
    # sh1 = flip(t(flip(t(sh0),1)),1)
    plot(make.raster(sh, d_flip), col = c("black", NA), add=T, legend=F)
  }
  # create final raster stack
  if (sw_totals){
    tfstk <- sw.totals.stk(keepAll = FALSE, mask=FALSE, savevar = NULL)
  } else{
    tfstk = sw.change.stk(zenith, average=TRUE, percentage=perc, savevar = NULL, mask_inner = FALSE)
  }
  print(proc.time() - ptm)
  return(tfstk)
}

###############################################
sw.moment <- function(sv, zenith, azimuth_eq, jd, alphaT = 0.45, plot_moment=NULL, zlen= NULL){
  # ! call topo.forcing.model
  # calculates all shortwave fluxes during the given moment of the day
  
  # zenith and incident angles
  cos_inc <- cos.slope(zenith, azimuth_eq, aspect = s_a[[2]], slope = s_a[[1]])
  cos_sfc <- cos(radians(zenith))
  
  # insolation arriving perpendicular to solar beam (direct and diffuse)
  Idirdif = insolation(zenith,jd,height,visibility,RH,tempK,0.002,0.45)
  Ib = matrix(Idirdif[,1],nrow=nrow(dem),ncol=ncol(dem))
  Id = matrix(Idirdif[,2],nrow=nrow(dem),ncol=ncol(dem))
  
  # terrain-reflected
  # Ir = Iglob * (1 - VF_mat) * alphaT) 
  # Ir = (Ib*sh + Id * VF_mat) * 0.45 * (1 - VF_mat) * cos_sfc[m]
  Iglob = (Ib + Id * VF_mat)*cos_sfc
  Ir = Iglob * (as.matrix((1 + cos(s_a[[1]]))/2) - VF_mat) * alphaT # (Hetrick 1993)
  Ir[Ir < 0 ] = 0
  # Ir = Iglob * (1 - VF_mat) * alphaT # Dozier (Hetrick 1993)
  
  # topographic shading
  sh <<- doshade(d_mat, sv, dl=30)
  # sh <<- doshade(as.matrix(d_flip), sv, dl=30)
  
  # plot moment (optional)
  if (!is.null(plot_moment)){
    p.sw.moment(Ib, Id, Ir, sh, cos_inc, cos_sfc, moment = plot_moment, zz = zlen)
  } else{
    # Run toposol models (values continuously saved in memory)
    topo.forcing.models(Ib, Id, Ir, sh, cos_inc, cos_sfc)
  }
}

###############################################
topo.forcing.models <- function(Ib, Id, Ir, sh, cos_inc, cos_sfc){
  # add values for moment 
  
  ## MODELS
  # model_flat      <<- model_flat        +   (Ib)*cos_sfc
  # model_flat_base <<- model_flat_base   +   (Ib + Id)*cos_sfc
  # model_flat_sh   <<- model_flat_sh     +   (Ib * sh)*cos_sfc
  # model_inc_sr    <<- model_inc_sr      +   (Ib) * cos_inc
  # model_inc_sh    <<- model_inc_sh      +   (Ib * sh)*cos_inc
  model_vf_base   <<- model_vf_base     +   (Ib * sh + Id)*cos_inc
  model_vf        <<- model_vf          +   (Ib * sh + Id * VF_mat )*cos_inc
  model_ref       <<- model_ref         +   (Ib * sh + Id * VF_mat + Ir)*cos_inc 
  
  print(paste("Ib:",round(mean(Ib,na.rm=T),2),' - Id:',round(mean(Id,na.rm=T),2),
              ' - Ir:',round(mean(Ir,na.rm=T),2)))
}
