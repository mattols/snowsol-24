#
#
#
#



# READ in HRRR
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/hrrr.t13z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210328/hrrr.t07z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210328/hrrr.t10z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/skiles-group1/HRRR_water_years/2021/hrrr.20210328/hrrr.t10z.wrfsfcf06.grib2'
grib_path = '/uufs/chpc.utah.edu/common/home/skiles-group2/HRRR_GSLB/hrrr.20220401/hrrr.t14z.wrfsfcf06.grib2'

grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210401/hrrr.t09z.wrfsfcf06.grib2'
topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
read.hrrr(grib_path, dpoly, which.hrrr=1)
sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
clear.sky.sw(dem, dpoly, dtime, savepath=NULL) # returns sw3
hrrr.sw.topo(hsw, cld, dem, savepath=NULL) # hsw3
snotel.loc(irwin_path)


## topofiles
topo.bounds(dpath)
s_a = slope.aspect(dem)
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 135)

sv_name = "workflow-terrain"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
dev.off()

sv_name = "workflow-viewf"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
plot(viewf, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
dev.off()




sv_name = "workflow-1"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(mask(h[['dswrf']],dpoly), add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
plot(h[['dswrf']], legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text="", side=3, font=2, line=2.5, cex=1.2))
dev.off()

# legend("bottomleft", c("DSWRF (h)", dtime), bty='n')


plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(hsw, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
plot(hsw, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))
# legend("bottomleft", c("DSWRF (hsw)"), bty='n')

sv_name = "workflow-2"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(hsw_, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
plot(hsw_, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text='', side=3, font=2, line=2.5, cex=1.2))
dev.off()
# legend("bottomleft", c("DSWRF (hsw)"), bty='n')

sv_name = "workflow-final"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(hsw3, add = T, col = heat.colors(100, 0.5), axes=F, legend = F)
plot(hsw3, legend.only=TRUE, col=heat.colors(100, 0.5),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text='', side=3, font=2, line=2.5, cex=1.2))
dev.off()
# legend("bottomleft", c("DSWRF (hsw3)"), bty='n')



sv_name = "workflow-k"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
new_col = colorspace::heat_hcl(50, c=c(80,30), l=c(30,90), power=1)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(k, add = T, col = new_col, axes=F, legend = F)
plot(k, legend.only=TRUE, col=new_col,
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text='', side=3, font=2, line=2.5, cex=1.2))
dev.off()


sv_name = "workflow-direct"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
new_col = heat.colors(100, 0.5) #colorspace::heat_hcl(50, c=c(80,30), l=c(30,90), power=1)
hsw_p = hsw1
hsw_p[1] = 350
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(hsw_p, add = T, col = new_col, axes=F, legend = F)
plot(hsw_p, legend.only=TRUE, col=new_col,
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text='', side=3, font=2, line=2.5, cex=1.2))
dev.off()

sv_name = "workflow-diffuse"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
new_col = heat.colors(100, 0.4) #colorspace::heat_hcl(50, c=c(80,30), l=c(30,90), power=1)
hsw_p = hsw02
hsw_p[1] = 350
plot(hsh, col=grey(c(0:100)/100), legend=FALSE,axes = 0,frame.plot=0, box = FALSE)
plot(hsw_p, add = T, col = new_col, axes=F, legend = F)
plot(hsw_p, legend.only=TRUE, col=new_col,
     legend.width=2, legend.shrink=0.5,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text='', side=3, font=2, line=2.5, cex=1.2))
dev.off()








# Assuming hsh, h, hsw, and hsw3 are pre-loaded raster objects
# and dtime is a variable for the legend
dswrf <- resample(h[['dswrf']], dem, method='ngb') #
# Set up the plotting layout (3 columns, 1 row)
par(mfrow=c(1, 3), mar=c(0, 0, 0, 0), oma = c(0,3,0,3))  # No margin between the plots
raster_list <- list(h[['dswrf']], hsw, dswrf)
min_val <- min(cellStats(stack(raster_list), min), na.rm = TRUE)  # Minimum value across all rasters
max_val <- max(cellStats(stack(raster_list), max), na.rm = TRUE)  # Maximum value across all rasters

# Define the color scale and breaks
col_scale <- heat.colors(100, 0.5)  # Set consistent color scheme
breaks <- seq(min_val, max_val, length.out = 21)  # Create breaks based on the min/max range


# Define the color scale for all plots (you can change this range to match your data)
col_scale <- heat.colors(100, 0.5)  # Set consistent color scheme
# First plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(h[['dswrf']], add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, breaks=breaks)
# Second plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(dswrf, add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, breaks=breaks)
# Third plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(hsw, add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, breaks=breaks)
# Add the legend (only once on the right side)
plot(h[['dswrf']], legend.only=TRUE, col=col_scale, breaks=breaks, 
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=bquote(W~m^2), side=3, font=2, line=2.5, cex=1.2))

par(mfrow=c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))


####
# PLOT SIDE-BY-SIDE
#
ngb <- resample(h[['dswrf']], dem, method='ngb') #
bilinear <- resample(h[['dswrf']], dem, method='bilinear') #
# CALC min max
min_dswrf <- cellStats(ngb, min, na.rm = TRUE)
max_dswrf <- cellStats(ngb, max, na.rm = TRUE)
min_hsw <- cellStats(bilinear, min, na.rm = TRUE)
max_hsw <- cellStats(bilinear, max, na.rm = TRUE)
min_hsw3 <- cellStats(hsw_, min, na.rm = TRUE)
max_hsw3 <- cellStats(hsw_, max, na.rm = TRUE)
# Get the global min and max values across all rasters
global_min <- min(min_dswrf, min_hsw, min_hsw3)
global_max <- max(max_dswrf, max_hsw, max_hsw3)



### NEW FIGURE
sv_name = "resample-methods"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=10, height=5, units = "in", res=300)

# Define the color scale and breaks
# Set up the plotting layout (3 columns, 1 row)
par(mfrow=c(1, 3), mar=c(0, 0, 0, 3), oma = c(0,5,0,5))  # No margin between the plots
col_scale <- heat.colors(100, 0.5)  # Set consistent color scheme
# First plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(ngb, add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, zlim=global_min, zmax=global_max)
# Second plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(bilinear, add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, zlim=global_min, zmax=global_max)
# Third plot
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(hsw_, add=TRUE, col=col_scale, axes=FALSE, legend=FALSE, zlim=global_min, zmax=global_max)
# Add the legend (only once on the right side)
plot(hsw_, legend.only=TRUE, col=col_scale, zlim=global_min, global_max, 
     legend.width=3.5, legend.shrink=0.50,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text="",side=3, font=2, line=2.5, cex=1.4))

dev.off()



par(mfrow=c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))



#######################################################
######### SAVE 




sv_name = "resample-compare"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
# DIFFERENCE from elevational resample
new_col = colorspace::diverge_hcl(100, c=100, l=c(50,90), power=1)
plot(hsw_ - bilinear, axes=F, box=FALSE, legend = F, 
     col= new_col)
plot(hsw_ - bilinear, legend.only=TRUE, col=new_col,
     legend.width=2.0, legend.shrink=0.50,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text = "",side=3, font=2, line=2.5, cex=1.4))
dev.off()

sv_name = "new-model-compare"
png(file.path("/uufs/chpc.utah.edu/common/home/u1037042/Documents/isnobal_data/erw/figs", paste0(sv_name, ".png") ), 
    width=7, height=12, units = "in", res=300)
# DIFFERENCE in model
model_dif = hsw03 - hsw_
model_dif[model_dif < -300] = -300
# new_col = colorspace::diverge_hcl(30, c=100, l=c(50,90), power=1)
new_col = colorspace::sequential_hcl(255)
plot(model_dif, axes=F, box=FALSE, legend = F, 
     col= new_col)
plot(model_dif, legend.only=TRUE, col=new_col,
     legend.width=2, legend.shrink=0.50,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text = "",side=3, font=2, line=2.5, cex=1.4))
dev.off()


## SHOW WORKFLOW
# 1.  HRRR
# HSW
# 2. Calculate K
# 3. Show Direct and Diffuse Beam


new_col = heat.colors(100, 0.7) #colorspace::heat_hcl(100, 0.5)
plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(hsw3, axes=F, box=FALSE, legend = F, 
     col= new_col)
plot(hsw3, legend.only=TRUE, col=new_col,
     legend.width=2, legend.shrink=0.50,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text = "",side=3, font=2, line=2.5, cex=1.4))





# DIFFERENCE in model
new_col = colorspace::diverge_hcl(30, c=100, l=c(50,90), power=1)
plot(hsw03 - hsw3, axes=F, box=FALSE, legend = F, 
     col= new_col)
plot(hsw03 - hsw3, legend.only=TRUE, col=new_col,
     legend.width=3.5, legend.shrink=0.60,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text = "",side=3, font=2, line=2.5, cex=1.4))

mod3_sum = calc(hsw3, "sum")
hsw_sum = calc(hsw_, "sum")


new_col = colorspace::diverge_hcl(30, c=100, l=c(50,90), power=1)
plot(mod3_sum - hsw_sum, axes=F, box=FALSE, legend = F, 
     col= new_col)
plot(mod3_sum - hsw_sum, legend.only=TRUE, col=new_col,
     legend.width=3.5, legend.shrink=0.60,
     axis.args=list(cex.axis=1.4),
     legend.args=list(text = "",side=3, font=2, line=2.5, cex=1.4))


# CLOUDS!

hrrr.moment.wrapper(grib_path, dpath, which.hrrr = which.hrrr, save_path = save_path)

par(mfrow=c(1,3))
topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
for (i in 8:19){
        if (i<10){
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210325/hrrr.t0',i,'z.wrfsfcf06.grib2')
        }else{
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210325/hrrr.t',i,'z.wrfsfcf06.grib2')
        }

        
        read.hrrr(grib_path, dpoly, which.hrrr=1) # returns h
        sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
        print(dtime)
        clear.sky.sw(dem, dpoly, dtime, savepath=NULL) # returns sw3
        hrrr.sw.topo(hsw, cld, dem, savepath=NULL) # hsw3, dswrf
        snotel.loc(irwin_path)
        
        # plot(k, main = dtime, col=blues9)
        # plot(new_k, main = zenith, col=blues9)
        # plot(h[['cld']], col=blues9, main =i)
        bilinear <- resample(h[['dswrf']], dem, method='bilinear') #
        min_dswrf <- cellStats(bilinear, min, na.rm = TRUE)
        max_dswrf <- cellStats(bilinear, max, na.rm = TRUE)
        min_hsw <- cellStats(hsw, min, na.rm = TRUE)
        max_hsw <- cellStats(hsw, max, na.rm = TRUE)
        global_min <- min(min_dswrf, min_hsw)
        global_max <- max(max_dswrf, max_hsw)
        
        col_scale <- heat.colors(100, 0.5) 
        new_col = colorspace::diverge_hcl(30, c=100, l=c(50,90), power=1)
        plot(bilinear, main = "bilinear", col=col_scale, zlim=c(global_min, global_max) )
        plot(hsw, main = "elevational", col=col_scale, zlim=c(global_min, global_max))
        plot(bilinear - hsw, col=new_col, main = "bilinear - hsw")
        
        # plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE, main = i)
        # plot(dswrf* cos_inc, col=heat.colors(255, 0.5), add=T);plot(pts, add=T,pch=1)
        # plot(pts, add=T,pch=1)
        # plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE, main = i)
        # plot(hsw3, col=heat.colors(255, 0.5), add=T);plot(pts, add=T,pch=1)
        # plot(pts, add=T,pch=1)
        # 
        # plot(values(dswrf)~values(sw0) )
        # points(values(dswrf*cos_inc), values(sw3), col='green')
        
}
par(mfrow=c(1,1))


par(mfrow=c(1,3))
mmdd = '0325'
for (i in 8:19){
        if (i<10){
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.2021',mmdd, '/hrrr.t0',i,'z.wrfsfcf06.grib2')
        }else{
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.2021',mmdd,'/hrrr.t',i,'z.wrfsfcf06.grib2')
        }
        
        read.hrrr(grib_path, dpoly, which.hrrr=1)
        #
        dat2 <- data.frame(x = values(h[['gmp']]), y =values(h[['dswrf']])) %>% na.omit()
        model <- lm(y ~ x, dat2)
        intercept <- coef(model)[1]  # Intercept (b0)
        slope <- coef(model)[2]     # Slope (b1)
        slope <- pmax(pmin(slope, 1), 0)
        intercept <- pmax(pmin(intercept, 1000), 0)
        #
        hd <- raster(dem); hd[] <- 1
        ff <-  res(h)[1]/res(dem)[1]
        hdis  <-  disaggregate(h[['dswrf']], fact = ff) # 3000/50 = 60
        hdd <-  resample(hdis, hd, method='ngb') # will not need to interpolate much but matches grid
        hsw <<- intercept + slope * dem
        
        # Correct for outlying points that could be cloud
        residuals_ = residuals(model)
        sigma_hat =  sd(residuals_) # check for other outliers?
        # Identify outliers: residuals > 3 standard deviations away from 0
        outlier_threshold <- 3 * sigma_hat
        outliers <- which(abs(residuals_) > outlier_threshold)
        # (optional) Cook's distance
        cooks_dist <- cooks.distance(model)
        influential_points <- which(cooks_dist > 4/length(dat2$x))
        #
        
        print(paste("for", i, "is clear?", slope > 0.00001))
        print(paste("R2", summary(model)$r.square))
        print(paste("pval", summary(model)$coefficients[2, 4] ))
        
        plot(h[['dswrf']], col=heat.colors(100, 0.5) , main =i)
        plot(intercept + slope * dem, col=heat.colors(100, 0.5), main = dtime)
        
        # line plot
        plot(y~x, data=dat2, main = dtime) 
        abline(model)
        # Plot points with outliers highlighted
        points(dat2$x[outliers], dat2$y[outliers], col='red')
        points(dat2$x[influential_points], dat2$y[influential_points], col='green')
        
        
}




par(mfrow=c(1,2))
topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
mmdd = '0325'
mmdd = '0401'
for (i in 8:19){
        if (i<10){
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.2021',mmdd, '/hrrr.t0',i,'z.wrfsfcf06.grib2')
        }else{
                grib_path = paste0('/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.2021',mmdd,'/hrrr.t',i,'z.wrfsfcf06.grib2')
        }
        
        read.hrrr(grib_path, dpoly, which.hrrr=1)
        
        
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
        if(slope > 0){
                # CALCULATE ELEVATIONAL DEPENDENCY
                hd <- raster(dem); hd[] <- 1
                ff <-  res(h)[1]/res(dem)[1]
                hdis  <-  disaggregate(h[['dswrf']], fact = ff) # 3000/50 = 60
                hdd <-  resample(hdis, hd, method='ngb') # will not need to interpolate much but matches grid
                hsw <<- intercept + slope * dem
                
                # Correct for outlying points that could be cloud
                residuals_ = residuals(model)
                sigma_hat =  sd(residuals_) # check for other outliers?
                # Identify outliers: residuals > 3 standard deviations away from 0
                outlier_threshold <- 3 * sigma_hat
                outliers <- which(abs(residuals_) > outlier_threshold)
                # (optional) Cook's distance
                cooks_dist <- cooks.distance(model)
                influential_points <- which(cooks_dist > 4/length(dat2$x))
                # This threshold for Cook's distance (4/n) is a common rule of thumb for identifying points that may be disproportionately influencing the model
                #
                # identify pixels
                hmk = raster(h) ;hmk[] <- 1;hmk[is.na(h[[1]])] <- NA
                hmk[as.numeric(names(outliers))] <- h[['dswrf']][as.numeric(names(outliers))] 
                hkp = rasterToPolygons(hmk)
                plot(hkp, col='red')
                
                hmk[as.numeric(names(outliers))] <- 0 
                hmkdis  <-  disaggregate(hmk, fact = ff)
                hmk_r  <-  resample(hmk, hdd)
                hmk0 <-  resample(hmkdis, hd, method='bilinear')
                plot(hmkdis)
                plot(hmk0)
                
                
                hmk = h[['dswrf']]
                hmk[as.numeric(names(outliers))] <- 
                hmkdis  <-  disaggregate(hmk, fact = ff)
                
                hdd2 = hdd
                plot(hdd2*hmkdis)
                hdd[match(values(hdd), h[['dswrf']][as.numeric(names(outliers))]) ]
                
                r1_modified <- overlay(hsw, hmk0, fun = function(x, y) { ifelse(!is.na(y), y, x) })
                
                hsw2 = hsw
                hsw2[which(!is.na(values(hmk)))] = hmk[!is.na(values(hmk))]
                plot(hsw2)
                
                
        }else{
                ### ELSE USE 
                hsw <<- resample(hdis, hd, method='bilinear')
        }
        
        
        # standard error
        sm$coefficients[1,2]
        
        plot(h[['dswrf']], col=heat.colors(100, 0.5) , main =i)
        plot(hmk)
        
        # line plot
        plot(y~x, data=dat2, main = dtime) 
        abline(model)
        # Plot points with outliers highlighted
        points(dat2$x[outliers], dat2$y[outliers], col='red')
        points(dat2$x[influential_points], dat2$y[influential_points], col='green')
        
        
        
        # bilinear <- resample(h[['dswrf']], dem, method='bilinear') #
        # min_dswrf <- cellStats(bilinear, min, na.rm = TRUE)
        # max_dswrf <- cellStats(bilinear, max, na.rm = TRUE)
        # min_hsw <- cellStats(hsw, min, na.rm = TRUE)
        # max_hsw <- cellStats(hsw, max, na.rm = TRUE)
        # global_min <- min(min_dswrf, min_hsw)
        # global_max <- max(max_dswrf, max_hsw)
        # 
        # col_scale <- heat.colors(100, 0.5) 
        # new_col = colorspace::diverge_hcl(30, c=100, l=c(50,90), power=1)
        # plot(bilinear, main = "bilinear", col=col_scale, zlim=c(global_min, global_max) )
        # plot(hsw, main = "elevational", col=col_scale, zlim=c(global_min, global_max))
        # plot(bilinear - hsw, col=new_col, main = "bilinear - hsw")
        

}






grib_path = '/uufs/chpc.utah.edu/common/home/uvu-group1/olson/snow-data/HRRR/hrrr.20210325/hrrr.t11z.wrfsfcf06.grib2'
topo.bounds(dpath) # returns: dem, viewf, slp, dpoly
sw.grid.elevation.resample(h, dem, method='bilinear') # returns hsw
clear.sky.sw(dem, dpoly, dtime, savepath=NULL) # returns sw3
hrrr.sw.topo(hsw, cld, dem, savepath=NULL) # hsw3, DSWRF
snotel.loc(irwin_path)

plot(hsh, col=grey(c(0:100)/100), legend=FALSE, axes=FALSE, frame.plot=FALSE, box=FALSE)
plot(dswrf, col=heat.colors(255, 0.5), add=T);plot(pts, add=T,pch=1)


plot(dswrf* cos_inc, col=heat.colors(255), main=names(dswrf_stk[[4]]))

