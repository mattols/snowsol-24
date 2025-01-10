#
# Snow depth from cast shadows
# test 01
#


ldem = raster("~/Documents/timp-ls-data/UTLC/snex/SNEX20_QSI_DEM_3M_USUTLC_20210921_20210921.tif")
lsd = raster("~/Documents/timp-ls-data/UTLC/snex/SNEX20_QSI_SD_3M_USUTLC_20210318_20210318.tif")

plot(ldem)
ldem


# unzip
# unzip("~/Documents/timp-ls-data/UTLC/sentinel/20210308.zip",
#       exdir = "~/Documents/timp-ls-data/UTLC/sentinel/20210308/")
# unzip("~/Documents/timp-ls-data/UTLC/sentinel/20210328.zip",
#       exdir = "~/Documents/timp-ls-data/UTLC/sentinel/20210328/")
# unzip("~/Documents/timp-ls-data/UTLC/planet/UTLC_202103_psscene_analytic_udm2.zip",
#       exdir = "~/Documents/timp-ls-data/UTLC/planet/202103/")
#
#
# # # # # # # #
# planet_scene
pth0 = "~/Documents/timp-ls-data/UTLC/planet/202103/PSScene/"
list.files(pth0,pattern="*3B_AnalyticMS_clip.tif") # raw data
list.files(pth0,pattern="*udm2_clip.tif") # classified?
list.files(pth0,pattern="3B_AnalyticMS_metadata_clip.xml")

# # # #
# create mosaic of 03-19-2021
ppath = list.files(pth0,pattern="20210319_.*3B_AnalyticMS_clip.tif",
           full.names = TRUE)
im19 = list(stack(ppath[1]), stack(ppath[2]))
im19 = do.call(raster::merge, im19)
im19
plotRGB(im19, 3,2,1)

# area poly
# lx = !is.na(ldemp)
# lx[lx==0]=NA
# dext = raster::rasterToPolygons(lx, n=8, dissolve=T)

# roi1
roi1 = click.new.roi(2e3,proj4string(im19))
im_crop = crop(im19, roi1)
plotRGB(im_crop, 3,2,1, stretch="lin")
ldemp = projectRaster(ldem, crs=raster::crs(proj4string(im19)))
ld = crop(ldemp, roi1)
plot(ld)

# snow
sd0 = projectRaster(lsd, crs=raster::crs(proj4string(im19)))
sd1 = crop(sd0, roi1)
plot(sd1)

# dem + sd
ld2 = ld + sd1

# meta
meta <- list.files(pth0,pattern="20210319_.*3B_AnalyticMS_metadata_clip.xml",
           full.names = TRUE)[1] %>% 
  xml2::read_xml() 
# eop:acquisitionDate
# ps:Acquisition

meta %>% xml2::xml_find_all("//ps:Acquisition")%>%xml2::xml_text()
sol_az = meta %>% xml2::xml_find_first("//opt:illuminationAzimuthAngle") %>%
  xml2::xml_text() %>% as.numeric()
sol_elv = meta %>% xml2::xml_find_first("//opt:illuminationElevationAngle") %>%
  xml2::xml_text() %>% as.numeric()
sol_vect <- insol::normalvector(90-sol_elv,sol_az)
# sol_vect2 <-  c(-sol_vect[1:2],sol_vect[3]))

shd <- insol::doshade(ld,sol_vect)
# sh_ls <- flip(t(flip(t(sh_ls),1)),1)

# test
plot(ld)
plot(shd, add=T, col=c("green",NA))
plotRGB(im_crop, 3,2,1)
plot(shd, add=T, col=c("green",NA))

# edit shade
# shd2 = shd*!is.na(ld)
shd2 = shd;shd2[is.na(ld)] = NA
# shd2 = !shd2; shd2[shd2==0] = NA
plot(shd2)
#
# filter pixels based on 'shade_size_filter'
# planet
shade_size_filter = 500; window_size=25
shclump <-  clump(!shd2, directions=8)
f<-as.data.frame(freq(shclump))
exludeShade <- f$value[which(f$count <= shade_size_filter)]
shfilter <- shclump
shfilter[shclump %in% exludeShade] <- NA 
shf <- (!is.na(shfilter))*!is.na(ld)
plot(shf)
# polygonize
shf2 = shf;shf2[shf2==0] = NA
shp = raster::rasterToPolygons(shf2, n=8, dissolve=T)
# saveRDS(shp, "~/Documents/timp-ls-data/UTLC/results/UTLC_shdpoly_500_20210318.rds")
# saveRDS(roi1, "~/Documents/timp-ls-data/UTLC/results/UTLC_roi1.rds")
# saveRDS(roi2, "~/Documents/timp-ls-data/UTLC/results/UTLC_roi2.rds")

# saveRDS(shp2, "~/Documents/timp-ls-data/UTLC/results/UTLC_shdpoly_plusSnow_20210318.rds")
# writeRaster(shf3, "~/Documents/timp-ls-data/UTLC/results/UTLC_SHADE_plusSnow_20210318.grd")

# quick save (Full image)
# writeRaster(shf2, "~/Documents/timp-ls-data/UTLC/results/UTLC_SHADE_500_20210318.grd")

# final plot
plotRGB(im_crop, 3,2,1)
plot(shp, add=T, border = "red")

# closer zoom
roi2 = click.new.roi(200,proj4string(im_crop))
imc2 = crop(im_crop, roi2)
plotRGB(imc2, 3,2,1)
plot(shp, add=T, border = "red")

# 3 meter res
# detect green threshold
# could use
plotRGB(imc2, 3,2,1)
spectra5 = raster::click(imc2, n = 5, cell=T)
names(spectra5) = c("cn", "blue", "green", "red", "nir")
spectra5 %>% reshape(varying=2:5,
                     direction = "long", v.names = "refl") %>% 
  ggplot2::ggplot(ggplot2::aes(x = time, y = refl, group = cn, colour=as.factor(id))) + ggplot2::geom_line() 

plotRGB(imc2, 3,2,1)
plot(imc2[[2]] < 14500, add=T, col=ggplot2::alpha(c("transparent","blue"), 0.1))
plot(shp, add=T, border = "red")

# snow depth
lsnow = projectRaster(lsd, crs=raster::crs(proj4string(im19)))
lsn = crop(lsnow, roi2)
plot(lsn)
plot(shp, add=T, border = "red")

plot(getValues(lsn))
# calculate distance 


# extract visual shadow
plotRGB(imc2)
plot(imc2[[2]]<max(imc2[[2]])*0.1)
perc10 = imc2[[2]]@data@max * 0.1
plot(imc2[[2]] < perc10)
# writeRaster(imc2, "~/Documents/timp-ls-data/UTLC/results/UTLC_Planet_Roi2_20210318.grd")

# entire one
roi3 = click.new.roi(400,proj4string(im_crop))
shx = crop(shf2, roi3)
shp3 = crop(shp, roi3)
plot(shx)
plot(shp3)
#


# HULL
shp_hull = gConvexHull(shp3)
plot(shp_hull)

# get rid of holes
outerRings = Filter(function(f){f@ringDir==1},shp3@polygons[[1]]@Polygons)
outerBounds = SpatialPolygons(list(Polygons(outerRings,ID=1)))
lines(outerBounds, col='green')
lines(shp3, col='red', lty=3)


# rotation
# rotang = 45
rotang = (sol_az - 180) # for mornings
rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center
center <- st_centroid(st_union(inpoly))

poly = st_as_sf(shp3) # %>% st_transform(4326)
inpoly = poly[1,] %>% sf::st_geometry()

tran(inpoly, rotang, center) %>% plot()
inpoly %>% plot(border="green",add=T)

oo = tran(inpoly, rotang, center)
#Extract coordinates
# bbox_list <- lapply(st_geometry(oo), st_bbox)
bbox = st_bbox(oo)

#To df
maxmin <- as.data.frame(matrix(unlist(bbox_list),nrow=nrow(oo)))
#get names
names(maxmin) <- names(bbox_list[[1]])
#Final shape
myshape2=bind_cols(oo,maxmin)
divline 
st_split(oo,divline) 

# x = 439964.1; y = 4487479


# # # # # # # #
## OLD OLD OLD
# # create mosaic from clip tile list
# r.mosaic <- lapply(list.files(pth0,pattern="*8b_clip.tif", full.names = TRUE), 
#                    terra::rast) %>% 
#   sprc %>% # create spatial raster list
#   mosaic
# r.mosaic
# # plot
# plotRGB(r.mosaic, 6,4,2, stretch='lin')
# # terra::writeRaster(r.mosaic, "D:/new23_data/planet/timp/planet/Mosaics/Planet8b_SR_20220918.tif")
# 














### BEAERING LINE
# my reproducible example
# https://stackoverflow.com/questions/64714899/r-find-distance-between-line-and-point-at-fixed-bearing-angles
library(sf)

# matrix of lon lat for the definition of the linestring
m<-rbind(
  c(12.09136, 45.86471),
  c(12.09120, 45.86495),
  c(12.09136, 45.86531),
  c(12.09137, 45.86540),
  c(12.09188, 45.86585),
  c(12.09200, 45.86592),
  c(12.09264, 45.86622),
  c(12.09329, 45.86624),
  c(12.09393, 45.86597),
  c(12.09410, 45.86585),
  c(12.09423, 45.86540),
  c(12.09411, 45.86495),
  c(12.09393, 45.86471),
  c(12.09383, 45.86451),
  c(12.09329, 45.86414),
  c(12.09264, 45.86413),
  c(12.09200, 45.86425),
  c(12.09151, 45.86451),
  c(12.09136, 45.86471)
)

# define the linestring
ls<-st_linestring(m)

# create a simple feature linestring with appropriate crs
ls<-st_sfc(ls, crs=4326)

# and now again going through the very same 
# definition process for a point

# define the origin point 
pt <- st_point(c(12.09286,45.86557))

# create simple feature point with appropriate crs
pt<-st_sfc(pt, crs = 4326)

plot(ls)
plot(pt, add=TRUE)

# get minimum distance from the origin point to the line
dist_min<-st_distance(ls, pt)

# get cordinates of the origin point
pt_orig<-st_coordinates(pt)

# load library for later use of the function destPoint()
library(geosphere)

# create vector of bearing angles of 10 degress amplitude
b_angles<-seq(0, 350, 10) 

# create empty container for final result as data frame
result<-data.frame(bearing=NULL, distance=NULL)

for(i in 1:length(b_angles)){
  
  result[i,"bearing"]<-b_angles[i]
  
  # calculate destination point coordinates with bearing angle i
  # at fixed safe distance (i.e. 100 times the minimum distance)
  # so that to avoid null intersection in next step calculation
  pt_dest<-destPoint(p=pt_orig, b=b_angles[i],d=dist_min*100)
  
  # define linestring from origin to destination
  b_ls<-st_sfc(st_linestring(rbind(pt_orig, pt_dest)), crs=4326)
  
  # get the intersection point between two features
  pt_int<-st_intersection(ls, b_ls)
  
  # get the distance
  d<-st_distance(pt, pt_int)
  
  result[i,"distance"]<-d
}




### PLOT BEARING
# https://stackoverflow.com/questions/22977453/plot-distance-and-bearing-in-r

# Define function to calculate coordinates given distance and bearing
get.coords <- function(a, d, x0, y0) {
  a <- ifelse(a <= 90, 90 - a, 450 - a)
  data.frame(x = x0 + d * cos(a / 180 * pi), 
             y = y0+ d * sin(a / 180 * pi))
}

# Set up plotting device
plot.new()
par(mar=c(2, 0, 0, 0), oma=rep(0, 4))
plot.window(xlim=c(-1100, 1100), ylim=c(-100, 1100), asp=1)

# Semicircles with radii = 100 through 1000
sapply(seq(100, 1000, 100), function(x) {
  lines(get.coords(seq(270, 450, length.out=1000), x, 0, 0))
})

# Horizontal line
segments(-1000, 0, 1000, 0)

# 45-degree lines
apply(get.coords(c(360-45, 45), 1000, 0, 0), 1, 
      function(x) lines(rbind(x, c(0, 0)), lwd=2))

# Plot white curves over black curves and add text
sapply(seq(100, 1000, 100), function(x) {
  txt <- paste0(x, 'm')
  w <- strwidth(txt, cex=0.9)/2
  a <- atan(w/x)/pi*180
  lines(get.coords(seq(-a, a, length=100), x, 0, 0), 
        lwd=2.5, col='white')
  text(0, x, txt, cex=0.8)
})

# Add points
points(with(dat, get.coords(-bear, dist, 0, 0)), pch=20)

# Add triangle
polygon(c(0, -30, 30), c(-5, -55, -55), col='black')






# https://stackoverflow.com/questions/51282724/creating-a-regular-polygon-grid-over-a-spatial-extent-rotated-by-a-given-angle
library(sf)
rotang = 45
rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center
inpoly <- st_read(system.file("shape/nc.shp", package="sf"))[1,] %>% 
  sf::st_transform(3857) %>% 
  sf::st_geometry()
center <- st_centroid(st_union(inpoly))
grd <- sf::st_make_grid(tran(inpoly, -rotang, center), cellsize = 3000)
grd_rot <- tran(grd, rotang, center)
plot(inpoly, col = "blue")
plot(grd_rot, add = TRUE)



# https://stackoverflow.com/questions/12663263/dissolve-holes-in-polygon-in-r
# 







poly = st_as_sf(shp3) # %>% st_transform(4326)
inpoly = poly[1,] %>% sf::st_geometry()
center <- st_centroid(st_union(inpoly))
grd <- sf::st_make_grid(tran(inpoly, -rotang, center), cellsize = 3000)
grd_rot <- tran(grd, rotang, center)
plot(inpoly, col = "blue")
plot(grd_rot, add = TRUE)

tran(inpoly, -45, center) %>% plot()

# HULL
shp_hull = gConvexHull(shp3)
plot(shp_hull)

# get rid of holes
outerRings = Filter(function(f){f@ringDir==1},shp3@polygons[[1]]@Polygons)
outerBounds = SpatialPolygons(list(Polygons(outerRings,ID=1)))
lines(outerBounds, col='green')
lines(shp3, col='red', lty=3)












#
#
# shade height
#

sc1_path = list.files("~/Documents/timp-ls-data/LC08_L1TP_038032_20230925_20231002_02_T1/", 
                      full.names = T, pattern = "B[1-7].TIF")
scene1 = stack(sc1_path)
plotRGB(scene1, 4,3,2, stretch="lin")

new1 = click.new.roi(2e3,proj4string(scene1))
crop(scene1, new1) %>% plotRGB( 4,3,2, stretch="lin")

save_base_name <- file.path(dirname(dirname(sc1_path)[1]), "dem", "scene1")
saveRDS(new1, paste0(save_base_name,"_glacierPolygons.rds"))

mini.shade.list <- function(folder, save_base_name){
  # get DEM
  d = raster::crop(d,shape) %>%
    raster::projectRaster(.,crs=raster::crs(proj4string(scene1)))
  
  # d = raster::crop(d,shape) %>%
  #   raster::projectRaster(.,crs=raster::crs(proj4string(scene1)), filename = paste0(save_base_name,"_dem",type,"30.grd"))
  # 
  # shade
  d_flip <- flip(t(raster::flip(t(d),1)),1)
  # solar geom
  xm_ls <- lapply(grep("MTL.xml",list.files(tile_list, full.names = T),value=TRUE), 
                  function(m) m %>% xml2::read_xml() )
  selv_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_ELEVATION")%>%xml2::xml_text()%>%as.numeric() )
  saz_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_AZIMUTH")%>%xml2::xml_text()%>%as.numeric() )
  sv_ls <- lapply(1:length(selv_ls),function(i) insol::normalvector(90-selv_ls[[i]][1],saz_ls[[i]][1]))
  sv_ls_cor <- lapply(sv_ls, function(sv) c(-sv[1:2],sv[3]))
  
  sh_ls <- insol::doshade(d_flip,sv_ls_cor[[5]])
  sh_ls <- flip(t(flip(t(sh_ls),1)),1)
  # extend to full tile
}

mini.shade <- function(scene_path, dem){
  scene_folder = dirname(scene_path)
  d_flip <- flip(t(raster::flip(t(dem),1)),1)
  # solar geom
  xm_ls <- lapply(grep("MTL.xml",list.files(scene_folder, full.names = T),value=TRUE), 
                  function(m) m %>% xml2::read_xml() )
  selv_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_ELEVATION")%>%xml2::xml_text()%>%as.numeric() )
  saz_ls <- lapply(xm_ls, function(x) x%>%xml2::xml_find_all("//SUN_AZIMUTH")%>%xml2::xml_text()%>%as.numeric() )
  sv_ls <- lapply(1:length(selv_ls),function(i) insol::normalvector(90-selv_ls[[i]][1],saz_ls[[i]][1]))
  sv_ls_cor <- lapply(sv_ls, function(sv) c(-sv[1:2],sv[3]))
  
  sh_ls <- insol::doshade(d_flip,sv_ls_cor[[5]])
  sh_ls <- flip(t(flip(t(sh_ls),1)),1)
  return(sh_ls)
}

sh1 = mini.shade(sc1_path, d)

# decent shadows
crop(scene1, new1) %>% plotRGB( 4,3,2, stretch="lin")
plot(sh1, add=T, col=c("green", NA))


d1 = raster("/uufs/chpc.utah.edu/common/home/u1037042/Documents/timp-ls-data/dem/dem_1.tif")
plot(d1)

sh2 = mini.shade(sc1_path, d1)

crop(scene1, new1) %>% plotRGB( 4,3,2, stretch="lin")
plot(sh2, add=T, col=c("green", NA))


#
sc2_path = list.files("~/Documents/timp-ls-data/LC09_L1TP_038032_20240530_20240531_02_T1/", 
                      full.names = T, pattern = "B[1-7].TIF")
scene2 = stack(sc2_path)
plotRGB(scene2, 4,3,2, stretch="lin")

#
sh3 = mini.shade(sc2_path, d1)

crop(scene2, new1) %>% plotRGB( 4,3,2, stretch="lin")
plot(sh3, add=T, col=c("green", NA))



fn2 = list.files(dirname(sc2_path)[1], full.names = T)
meta2 = grep("MTL.xml",fn2, value = T)
metaxl = meta2 %>% xml2::read_xml() 

metaxl %>% 
  xml2::xml_find_all("//SUN_ELEVATION")%>%xml2::xml_text()%>%as.numeric()
