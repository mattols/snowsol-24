# Pat Solar!

# get variables
file_path = '/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/OLI_snow_scenes/LC08_L2SP_034033_20190105_20200829_02_T1_MTL.txt'
# lns = read.delim(file_path) # gsub("([0-9]+)\\..*$", "\\1", ll)
con = file(file_path)
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


# all
list.files('/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/OLI_snow_scenes', 
           pattern = "*MTL.txt")

# get b1 extent and download NasaDEM for each

list.files('/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/OLI_snow_scenes', 
           pattern = "20150126.*[3-5].TIF", full.names = T) %>% brick() %>% 
  plotRGB(3,2,1)

all_b1 = list.files('/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/OLI_snow_scenes', 
           pattern = "*B1.TIF", full.names = T)
plot(raster(all_b1[1]))




dem = raster('/uufs/chpc.utah.edu/common/home/skiles-group2/pnaple/DEM/dem90_hf_30m.tif')
plot(dem)

### PROBLEMS DOWNLOADING
# NEED MULTIPLE DEMs to MOSAIC
# function for single DEM download
# input rpath | output DEM location
# convert bottomright to latlon
points1 = rbind(bbox(raster(all_b1[1]))[,1])
sputm <- SpatialPoints(points1, proj4string=raster::crs(raster(all_b1[1])))
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
seq_x_y = floor(bbox(spgeo)[,1])
# lnlt <- coordinates(spgeo) # bbox(raster(all_b1[1]))
# Generate tile names from measure
dem_base = "http://e4ftl01.cr.usgs.gov/MEASURES/NASADEM_HGT.001/2000.02.11/"
tile_names <- paste0("N",gsub(pattern=" ","W", paste(seq_x_y[2],abs(seq_x_y[1]))))
dem_list <- paste0(dem_base,"NASADEM_HGT_",tolower(tile_names),".zip")
# download to tmp
if(!exists("username")){username<<-readline(prompt = "Earthdata email:");pass<<-readline(prompt = "Earthdata password:")} 
cat("...retrieving DEM", tile_names, "at \n", dem_list, "\n")
temp <- tempfile(pattern = tile_names, fileext = ".zip")
httr::GET(dem_list, httr::authenticate(username,pass,type='basic'),
          httr::write_disk(temp))
out <- unzip(temp, exdir = tempdir()); file.remove(temp)
rl <- lapply(list.files(tempdir(),full.names = T,pattern = ".hgt"), raster)
cat("Reprojecting and saving SRTM DEM mosaic...")
do.call(merge, c(rl, tolerance = 1)) %>% crop(.,shape) %>% 
  mask(.,shape,filename = paste0(save_base_name,"_demMOSAIC.grd")) 
raster(paste0(save_base_name,"_demMOSAIC.grd")) %>%
  raster::projectRaster(.,crs=landsat_crs, filename = paste0(save_base_name,"_dem",type,"30.grd"))
unlink(paste0(save_base_name,"_demMOSAIC.grd"));unlink(paste0(save_base_name,"_demMOSAIC.gri"))
unlink(list.files(tempdir(), full.names = T,pattern="^[Nn]"))
cat("\n file",paste0("'dem",type,"30.grd'"), "created")


raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/trash/NASADEM_HGT_n39w108/n39w108.hgt") %>% 
  plot(add=T)


library(RCurl)
url <- paste0(dem_base,"NASADEM_HGT_",tolower(tile_names),".zip")
content <- getBinaryURL(url, userpwd = paste0(username,":",pass), 
                        ftp.use.epsv = FALSE)
writeBin(content, con = basename(url))


find.dem <- function(save_base_name, type=c("NASADEM", "SRTM")){
  # download SRTM DEM 30-meter for larger tile
  landsat_crs = raster::crs(readRDS(paste0(save_base_name,"_glacierPolygons.rds")))
  shape = as(raster::extent(readRDS(paste0(save_base_name,"_glacierPolygons.rds"))), 'SpatialPolygons')
  raster::crs(shape) <- raster::crs(readRDS(paste0(save_base_name,"_glacierPolygons.rds")))
  shape = sp::spTransform(shape,raster::crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # generate list (NEW)
  seqX <- seq(floor(min(shape@polygons[[1]]@Polygons[[1]]@coords[,1])),
              floor(max(shape@polygons[[1]]@Polygons[[1]]@coords[,1])),1)
  seqY <- seq(floor(min(shape@polygons[[1]]@Polygons[[1]]@coords[,2])),
              floor(max(shape@polygons[[1]]@Polygons[[1]]@coords[,2])),1)
  lstiles <- unlist(lapply(1:length(seqY),function(y) lapply(1:length(seqX),function(x) paste(seqY[y],seqX[x])) ))
  tile_names <- paste0("N",gsub(pattern=" ","E0", lstiles))
  #
  if(type=="SRTM"){
    dem_base = "http://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/"
    dem_list <- paste0(dem_base,tile_names,".SRTMGL1.hgt.zip")
  }else if(type=="NASADEM"){
    dem_base = "http://e4ftl01.cr.usgs.gov/MEASURES/NASADEM_HGT.001/2000.02.11/"
    dem_list <- paste0(dem_base,"NASADEM_HGT_",tolower(tile_names),".zip")
  }else{stop("DEM type not recognized!")}
  if(!exists("username")){username<<-readline(prompt = "Earthdata email:");pass<<-readline(prompt = "Earthdata password:")} 
  cat("\n")  # "matthew.olson@geog.utah.edu" | "Televator9" # could move this to beginning of code
  for (t in 1:length(tile_names)){
    print(paste("...retrieving", t,"of",length(tile_names), type,"files."),q=F)
    temp <- tempfile(pattern = tile_names[t], fileext = ".zip")
    httr::GET(dem_list[t],httr::authenticate(username,pass,type='basic'),httr::write_disk(temp))
    out <- unzip(temp, exdir = tempdir()); file.remove(temp)
  }
  rl <- lapply(list.files(tempdir(),full.names = T,pattern = ".hgt"), raster)
  cat("Reprojecting and saving SRTM DEM mosaic...")
  do.call(merge, c(rl, tolerance = 1)) %>% crop(.,shape) %>% 
    mask(.,shape,filename = paste0(save_base_name,"_demMOSAIC.grd")) 
  raster(paste0(save_base_name,"_demMOSAIC.grd")) %>%
    raster::projectRaster(.,crs=landsat_crs, filename = paste0(save_base_name,"_dem",type,"30.grd"))
  unlink(paste0(save_base_name,"_demMOSAIC.grd"));unlink(paste0(save_base_name,"_demMOSAIC.gri"))
  unlink(list.files(tempdir(), full.names = T,pattern="^[Nn]"))
  cat("\n file",paste0("'dem",type,"30.grd'"), "created")
}



# Calculate insol

###References:
#Iqbal, M. (2012). An introduction to solar radiation. Elsevier.

###Inputs:
#  Esc=1367     [Wm-2]             (Solar constant)
#  sza          [radians]          (zenith_angle) 
#  press        [mb]               (local barometric)
#  albedo       [double]           (surface/ground albedo)
#  ang_alpha    [dimensionless]    (Angstrom_exponent, also known as alpha)
#  ang_beta     [dimensionless]    (Angstrom turbidity coefficient)
#  ozone        [atm.cm]           (total columnar amount of ozone)
#  wv           [atm.cm]           (total columnar amount of water vapour)

###Outputs:
#  Ebn  [Wm-2]     (Direct normal irradiance)
#  Edh  [Wm-2]     (Diffuse horizontal irradiance)
#  Egh  [Wm-2]     (Global horizontal irradiance) 

###Notes:
#  Dayth is the day number ranging from 1 to 365. 

###Codes:
IrradianceIqbalc<-function(){
  #extraterrestrial irradiance
  Esc=1367
  totaldayofyear=366-ceiling(Year/4-trunc(Year/4)) 
  B=(Dayth-1)*2*pi/totaldayofyear
  Eext=Esc*(1.00011+0.034221*cos(B)+0.00128*sin(B)+0.000719*cos(2*B)+0.000077*sin(2*B))
  
  #air mass
  am=(cos(sza)+0.15*(3.885+(pi/2-sza)/pi*180)^-1.253)^-1
  ama=am*press/1013.25
  
  #the transmittance by Rayleigh scattering
  TR=exp(-0.0903*(ama)^0.84*(1+ama-ama^1.01)) 
  
  #the transmittance by ozone
  U3=ozone*am       
  alpho=(0.1611*U3*(1+139.48*U3)^-0.3035 - 0.002715*U3*(1.0+0.044*U3+0.0003*U3^2)^-1)
  TO=1-alpho
  
  #the transmittance by uniformly mixed gases
  TG=exp(-0.0127*ama^0.26)
  
  #the transmittance by water vapor
  PW=am*wv*(press/1013.25)^0.75
  alphw=2.4959*PW*((1+79.034*PW)^0.6828+6.385*PW)^-1
  TW=1-alphw
  
  #the transmittance by Aerosol
  kaa=ang_beta*(0.2758*0.38^-ang_alpha+0.35*0.5^-ang_alpha)     
  TA=exp(-kaa^0.873*(1+kaa-kaa^0.7088)*ama^0.9108)
  
  #direct beam irradiance
  Ebniqbalc=0.9751*Eext*TR*TO*TG*TW*TA
  
  #the transmittance of direct radiation due to aerosol absorptance
  w0=0.9  #as suggested by author from Bird Hulstrom
  TAA=1-(1-w0)*(1-ama+ama^1.06)*(1-TA)
  #the Rayleigh-scattered diffuse irradiance after the first pass through the atmosphere
  Edr = 0.79*Eext*cos(sza)*TO*TG*TW*TAA*0.5*(1-TR)/(1-ama+ama^1.02)
  
  #the aerosol scattered diffuse irradiance after the first pass through the atmosphere
  TAS=TA/TAA
  Fc=0.84   #as suggested by author
  Eda=0.79*Eext*cos(sza)*TO*TG*TW*TAA*Fc*(1-TAS)/(1-ama+ama^1.02)
  
  #global horizontal irradiance
  rg=albedo   #ground albedo
  ra=0.0685+(1-Fc)*(1-TAS)
  Eghiqbalc=(Ebniqbalc*cos(sza)+Edr+Eda)/(1-rg*ra)
  #diffuse horizontal irradiance
  Edhiqbalc=Eghiqbalc-Ebniqbalc*cos(sza)
  
  ###Quality control
  lower=0
  Ebniqbalc[Ebniqbalc<lower]=0
  Edhiqbalc[Edhiqbalc<lower]=0
  Eghiqbalc[Eghiqbalc<lower]=0
  return(list(Ebniqbalc, Edhiqbalc, Eghiqbalc))
}