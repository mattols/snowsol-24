#
# R instant isolation
# R 4.03 Geospatial Packages (CHPC)
# 11/2024
#
# # # # # #

rm(list = ls());gc()
source('~/Documents/src/snowtools/Pat-olirf/landsat-insol-1.R', echo=FALSE)

# # # # # # #
# Workflow  #
# start
strt = Sys.time()

# create mtl paths
list.scene.info()

# TEST SINGLE
mtl_path = mtl_files[1]
terrain.variables(mtl_path, dpath_main)
extract.solar.vars(mtl_path)
insol.params() 
# calculate solar
make.shade2(90-s_elv, s_az) # 2.5+ mins
Sys.time() - strt 
i0 = sw.moment.mtl(mtl_path, save_sw = FALSE)
i1 = sw.solar.noon(mtl_path, save_sw = FALSE)
# end time
Sys.time() - strt 

# check zoom (area 1)
# pts = click(n=1)
pts = cbind(392236, 4394872)
pts2 <- SpatialPoints(coords = pts, proj4string = CRS(proj4string(d)) )
plot(pts2, add=T,col='red', cex=3, pch=20)
pb = buffer(pts2, 1.6e4)
plot(extent(pb), add=T, col='red')

## FINAL PLOTS
# RGB

# Visualize INSOL with hillshade
sw2 = i1; d_name = dmoment_char_noon #solar noon
sw2 = i0; d_name = dmoment_char # landsat moment
# plot
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(sw2, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
plot(sw2, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
legend("bottomleft", c("Clear-sky Modeled Insolation",
                       paste("TilePath:", strsplit(basename(mtl_path), "_")[[1]][3], "| Date:",
                             d_name ) ) )






################################
 # save to insol folder
# START LOOP
source('~/Documents/src/snowtools/Pat-olirf/landsat-insol-1.R', echo=FALSE)
list.scene.info()
for (i in 49:length(mtl_files)){ #length(mtl_files)){
  # start time
  strt = Sys.time()
  # variables
  list.scene.info() # create mtl paths
  print(paste("Starting clear-sky insol", i, "of", length(mtl_files)))
  mtl_path = mtl_files[i] # current file
  terrain.variables(mtl_path, dpath_main) 
  extract.solar.vars(mtl_path)
  insol.params() # 13 sec to here
  # calculate shade
  make.shade2(90-s_elv, s_az) # 2.5+ mins
  sw.moment.mtl(mtl_path, save_sw = TRUE) # 3.5 mins combined (for small)
  sw.solar.noon(mtl_path, save_sw = TRUE)
  print(Sys.time() - strt) #END
  # clear cache
  rm(list = setdiff(ls(), c("i", "mtl_files")));gc()
  source('~/Documents/src/snowtools/Pat-olirf/landsat-insol-1.R', echo=FALSE)
}

###############
# # # # # # # #
# Check RESULTS
source('~/Documents/src/snowtools/Pat-olirf/landsat-insol-1.R', echo=FALSE)
list.scene.info()
r_mom <- list.files("/uufs/chpc.utah.edu/common/home/u1037042/olson/snow-data/drfs/oli-drfs/agu24/insol2", 
           full.names = T, pattern = "moment")
# check that they match
gsub("^dsw_(.*)_moment.*", '\\1', basename(r_mom)) == all_dates[1:length(r_mom)]
paste("Ended with", length(r_mom))
 
ndat = gsub(".*moment_(.*).tif", '\\1', basename(r_mom))
plot(raster(r_mom[1]), col = heat.colors(100, 0.5), main=ndat[13])
plot(raster(r_mom[14]), col = heat.colors(100, 0.5),  main=ndat[14])

# solar noon
rsol <- list.files("/uufs/chpc.utah.edu/common/home/u1037042/olson/snow-data/drfs/oli-drfs/agu24/insol2", 
                    full.names = T, pattern = "solnoon")
ndat = gsub(".*solnoon_(.*).tif", '\\1', basename(rsol))
plot(raster(rsol[1]), col = heat.colors(100, 0.5), main=ndat[1])
plot(raster(rsol[13]), col = heat.colors(100, 0.5), main=ndat[13])
plot(raster(rsol[14]), col = heat.colors(100, 0.5),  main=ndat[14])
     





# Share with Pat (u1369303)
# sudo chown username:u1369303 /uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/drfs/oli-drfs/agu24/insol
# sudo chmod 755 /uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/drfs/oli-drfs/agu24/insol
# 
