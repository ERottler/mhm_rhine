###

#Analyze mhm model results

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales, ncdf4,
               ncdf4.helpers, sp, raster, viridis, meltimr, POT)

# run_dir <- "D:/nrc_user/rottler/mhm_run/6935053/"
run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"

bas_dir <- "U:/rhine_fut/R/"

#load functions
source(paste0(bas_dir, "mhm_rhine/functs.R"))

sta_yea <- 1954
end_yea <- 2013

stopCluster(my_clust)

n_cores <- 25 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr))
registerDoParallel(my_clust)

#inp_vis----

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

dem <- raster(paste0(run_dir, "input/morph/dem.asc"))
plot(dem, col = viridis(100))

# basins <-  rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp", encoding = "UTF8")
# basin <- spTransform(basins[basins@data$Ort == "Rheinfelden, Messstation",], CRS = crswgs84)
basin_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/lobith_catch.shp")
basin <- spTransform(basin_raw, CRS = crswgs84)

#Load ncdf E-OBS gridded datasets
nc_temp_file <- paste0(run_dir, "input/meteo/tavg.nc")
nc_prec_file <- paste0(run_dir, "input/meteo/pre.nc")
nc_petr_file <- paste0(run_dir, "input/meteo_HS/pet.nc")

nc_temp <- nc_open(nc_temp_file)
nc_prec <- nc_open(nc_prec_file)
nc_petr <- nc_open(nc_petr_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_temp, varid = "lon2D")
lat <- ncdf4::ncvar_get(nc_temp, varid = "lat2D")
date <- as.Date(as.character(nc.get.time.series(nc_temp, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-01")
count_date <- length(date) - sta_date_ind

temps_cube <- ncvar_get(nc_temp, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lat), count_date), varid = "tavg")

precs_cube <- ncvar_get(nc_prec, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lat), count_date), varid = "pre")

evapo_cube <- ncvar_get(nc_petr, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lat), count_date), varid = "pet")

temps_mea <- apply(temps_cube, c(1,2), mea_na)
precs_mea <- apply(precs_cube, c(1,2), sum_na) / length(sta_yea:end_yea)
evapo_mea <- apply(evapo_cube, c(1,2), sum_na) / length(sta_yea:end_yea)

cols_spat_tem <- foreach(i = 1:length(c(temps_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(temps_mea),
          dat_ref = c(temps_mea),
          do_bicol = F,
          virid_dir = 1)
  
}
cols_spat_pre <- foreach(i = 1:length(c(precs_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(precs_mea),
          dat_ref = c(precs_mea),
          do_bicol = F)
  
}
cols_spat_eva <- foreach(i = 1:length(c(evapo_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(evapo_mea),
          dat_ref = c(evapo_mea),
          do_bicol = F)
  
}



pdf(paste0(bas_dir,"res_figs/met_map.pdf"), width = 16, height = 6)

layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4, rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())

par(family = "serif")

cex_pch <- 1.32

#Plot Temperature
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin, xlim = c(range(c(lon))), ylim = c(range(c(lat))))
points(c(lon), c(lat), pch = 15, col = cols_spat_tem, cex = cex_pch)
par(new = T)
plot(basin, add = T)
mtext("a) Temperature", side = 3, line = -1.0, cex = 1.7)


par(mar = c(2.0, 0.2, 5.0, 2.9))

my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = 1)))(200))
my_bre <- seq(range(temps_mea)[1], range(temps_mea)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(sc_doy_simu), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[°C]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot Precipitation
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin, xlim = c(range(c(lon))), ylim = c(range(c(lat))))
points(c(lon), c(lat), pch = 15, col = cols_spat_pre, cex = cex_pch)
par(new = T)
plot(basin, add = T)
mtext("b) Precipitation", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))

my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(precs_mea)[1], range(precs_mea)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(sc_doy_simu), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot Evapotrationspiration
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin, xlim = c(range(c(lon))), ylim = c(range(c(lat))))
points(c(lon), c(lat), pch = 15, col = cols_spat_eva, cex = cex_pch)
par(new = T)
plot(basin, add = T)
mtext("c) Evapotranspiration", side = 3, line = -1.0, cex = 1.7)


par(mar = c(2.0, 0.2, 5.0, 2.9))

my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(evapo_mea)[1], range(evapo_mea)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(sc_doy_simu), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

dev.off()



#Land use maps
lu_1990 <- raster(paste0(run_dir, "input/luse/lc_1990_3class.asc"))
lu_2000 <- raster(paste0(run_dir, "input/luse/lc_2000_3class.asc"))
lu_2006 <- raster(paste0(run_dir, "input/luse/lc_2006_3class.asc"))
lu_2012 <- raster(paste0(run_dir, "input/luse/lc_2012_3class.asc"))


pdf(paste0(bas_dir,"res_figs/lu_map.pdf"), width = 16, height = 6)

cols_lu <- c("green4", "grey25", "gold")

par(family = "serif")
par(mfrow = c(1, 4))
par(mar = c(0.5, 2, 2, 2))

plot(lu_1990, col = cols_lu, axes=F, legend = F, ylab = "", xlab = "", box = F)
legend("topright", c("forest", "impervious", "pervious"), pch = 15, col = cols_lu, bty = "n", cex = 1.4)
mtext("1950-1995", side = 3, line = 0.0, cex = 1.2)

plot(lu_2000, col = cols_lu, axes=F, legend = F, ylab = "", xlab = "", box = F)
mtext("1996-2003", side = 3, line = 0.0, cex = 1.2)

plot(lu_2006, col = cols_lu, axes=F, legend = F, ylab = "", xlab = "", box = F)
mtext("2004-2009", side = 3, line = 0.0, cex = 1.2)

plot(lu_2012, col = cols_lu, axes=F, legend = F, ylab = "", xlab = "", box = F)
mtext("2010-2016", side = 3, line = 0.0, cex = 1.2)

dev.off()


#LAI

lai_class <- raster(paste0(run_dir, "input/morph/LAI_class.asc"))

lai_clas_1 <- c("Evergreen-Needleleaf-Forest", "darkgreen")
lai_clas_2 <- c("Evergreen-Broadleaf-Forest", "forestgreen")
lai_clas_3 <- c("Deciduous-Needleleaf-Forest", "green4")
lai_clas_4 <- c("Deciduous-Broadleaf-Forest", "chartreuse4")
lai_clas_5 <- c("Mixed-Forests", "darkolivegreen4")
lai_clas_6 <- c("Closed-Shrublands", "moccasin")
lai_clas_7 <- c("Open-Shrublands", "moccasin")
lai_clas_8 <- c("Woody-Savannas", "white")
lai_clas_9 <- c("Savannas", "white")
lai_clas_10 <- c("Grasslands", "chartreuse")
lai_clas_11 <- c("Permanent-wetlands", "darkslategray4")
lai_clas_12 <- c("Croplands", "darkred")
lai_clas_13 <- c("Urban-and-Built-Up", "grey45")
lai_clas_14 <- c("cropland-natural-vegetation-mosaic", "red2")
lai_clas_15 <- c("Snow-and-Ice", "cadetblue1")
lai_clas_16 <- c("Barren-or-Sparsely-Vegetated", "salmon1")
lai_clas_17 <- c("Water", "darkblue")
lai_clas_18 <- c("Wooded-Tundra", "white")
lai_clas_19 <- c("Mixed-Tundra", "white")
lai_clas_20 <- c("Barren-Tundra", "white")
lai_clas_21 <- c("SEALED-AREA-75%", "black")
lai_clas_22 <- c("SEALED-AREA-50%", "grey25")
lai_clas_23 <- c("SEALED-AREA-25%", "grey35")

lai_clas_col <- rbind(lai_clas_1, lai_clas_2, lai_clas_3, lai_clas_4, lai_clas_5, lai_clas_6, lai_clas_7, lai_clas_8, lai_clas_9,
                      lai_clas_10, lai_clas_11, lai_clas_12, lai_clas_13, lai_clas_14, lai_clas_15, lai_clas_16, lai_clas_17,
                      lai_clas_18, lai_clas_19, lai_clas_20, lai_clas_21, lai_clas_22, lai_clas_23)

pdf(paste0(bas_dir, "res_figs/lai_map.pdf"), width = 8, height = 6)

par(family = "serif")
par(mfrow = c(1, 2))

par(mar = c(0, 0, 0, 0))
plot(lai_class, breaks = seq(0.5, 23.5, 1), col = lai_clas_col[, 2], axes=F, legend = F, ylab = "", xlab = "", box = F)
# mtext("LAI classes", side = 3, line = -0.5)
par(mar = c(0, 1, 1, 1))
plot(1:100, 1:100, type = "n", ylab = "", xlab = "", axes = F)
points(rep(5, 23), seq(5, 95, length.out = 23), pch = 19, col = lai_clas_col[, 2])
for(i in 1:nrow(lai_clas_col)){
  text(x = 10, y = seq(5, 95, length.out = 23)[i], lai_clas_col[i, 1], adj = 0)
}

dev.off()


#Soil

soil_bas_1 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_01.asc"))
soil_bas_2 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_02.asc"))
soil_bas_3 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_03.asc"))
soil_bas_4 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_04.asc"))
soil_bas_5 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_05.asc"))
soil_bas_6 <- raster(paste0(run_dir, "input/morph/soil_class_horizon_06.asc"))

pdf(paste0(bas_dir, "res_figs/soil_map.pdf"), width = 12, height = 8)

par(mfrow = c(2, 3))
par(family = "serif")
par(mar = c(0, 0, 3, 5))


plot(soil_bas_1, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("0 - 50 mm ", side = 3, line = 0.0)

plot(soil_bas_2, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("51 - 150 mm", side = 3, line = 0.0)

plot(soil_bas_3, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("151 - 300 mm", side = 3, line = 0.0)

plot(soil_bas_4, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("301 - 600 mm", side = 3, line = 0.0)

plot(soil_bas_5, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("600 - 1000 mm", side = 3, line = 0.0)

plot(soil_bas_6, col = viridis(100, direction = -1), zlim = c(0, 2000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("1001 - 2000 mm", side = 3, line = 0.0)

dev.off()


#Further morphological input

slo_bas <- raster(paste0(run_dir, "input/morph/slope.asc"))
asp_bas <- raster(paste0(run_dir, "input/morph/aspect.asc"))
geo_bas <- raster(paste0(run_dir, "input/morph/geology_class.asc"))
fdir_bas <- raster(paste0(run_dir, "input/morph/fdir.asc"))
facc_bas <- raster(paste0(run_dir, "input/morph/facc.asc"))


pdf(paste0(bas_dir, "res_figs/morph_maps.pdf"), width = 12, height = 8)

par(mfrow = c(2, 3))
par(family = "serif")
par(mar = c(0, 0, 3, 5))

plot(dem, col = viridis(100, direction = -1), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Elevation", side = 3, line = 0.0)

plot(slo_bas, col = viridis(100, direction = -1), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Slope", side = 3, line = 0.0)

plot(asp_bas, col = viridis(100, direction = -1), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Aspect", side = 3, line = 0.0)

plot(geo_bas, breaks = seq(0.5, 8.5, 1), col = viridis(8), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Geology", side = 3, line = 0.0)

plot(fdir_bas, col = viridis(100, direction = -1), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Flow direction", side = 3, line = 0.0)

plot(facc_bas, col = c("grey92", "blue2", "darkblue"), breaks = c(0, 10000, 1000000, 10000000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Flow accumulation", side = 3, line = 0.0)

dev.off()



#dis_vis----

dis_mhm <- read.table(paste0(run_dir, "output/daily_discharge.out"), header = T)
dis_mhm$date <- as.Date(strptime(paste0(dis_mhm$Day, "-", dis_mhm$Mon, "-", dis_mhm$Year), "%d-%m-%Y", tz="UTC"))

#Runoff seasonality

quants <- seq(0.01, 0.99, by = 0.01)

f_qvalu_obs <- function(quant_sel){dis_ana(disc = dis_mhm$Qobs_0006435060,
                                           date = dis_mhm$date,
                                           start_year = sta_yea,
                                           end_year = end_yea,
                                           break_day = 0,
                                           quant_in = quant_sel,
                                           do_moving_average = F,
                                           window_width = 30,
                                           method_analys = "quantile",
                                           method_quant = "empirical"
)}

qvalu_obs <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

f_qvalu_sim <- function(quant_sel){dis_ana(disc = dis_mhm$Qsim_0006435060,
                                           date = dis_mhm$date,
                                           start_year = sta_yea,
                                           end_year = end_yea,
                                           break_day = 0,
                                           quant_in = quant_sel,
                                           do_moving_average = F,
                                           window_width = 30,
                                           method_analys = "quantile",
                                           method_quant = "empirical"
)}

qvalu_sim <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif <- qvalu_obs - qvalu_sim


#Plot: Seasonality of runoff
pdf(paste0(bas_dir, "res_figs/runoff_qu.pdf"), width = 8, height = 6)

par(family = "serif")

layout(matrix(c(rep(1, 8), 2,
                rep(5, 8), 6,
                rep(3, 8), 4),
              3, 9, byrow = T), widths=c(), heights=c())

cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_obs, qvalu_sim)), alptempr::max_na(c(qvalu_obs, qvalu_sim)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_obs, cols = my_col, breaks = my_bre, header = "a) Qobs", lab_unit = "[m³/s]")


cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_sim, qvalu_obs)), alptempr::max_na(c(qvalu_sim, qvalu_obs)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_sim, cols = my_col, breaks = my_bre, header = "c) Qsim", lab_unit = "[m³/s]")


cols_max <- colorRampPalette(c(rep("grey98", 10), "lightgoldenrod2", "gold3", "goldenrod3", "orangered4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1, 1,2,3,4)], rep("lightcyan3", 1), rep("grey98", 10)))(100)
my_col <- c(cols_min, cols_max)
# my_bre <- seq(-max_na(abs(qvalu_dif)), max_na(abs(qvalu_dif)), length.out = length(my_col)+1)

my_bre <- c(-(lseq(0.01, max_na(abs(qvalu_dif)), length.out = length(my_col)/2)[(length(my_col)/2):1]),
            lseq(0.01, max_na(abs(qvalu_dif)), length.out = length(my_col)/2+1))

dis_image(data_plot = qvalu_dif, cols = my_col, breaks = my_bre, header = "b) Qobs - Qsim", 
          lab_unit = "[m³/s]", do_cont = F)

dev.off()



#Plot: Runoff time series
pdf(paste0(bas_dir, "res_figs/runoff_ts.pdf"), width = 8, height = 4)

par(mfrow = c(3, 1))
par(mar = c(2, 3, 0.5, 0.5))

vis_run(sta_day = "1954-01-01",
        end_day = "1956-12-31",
        do_legend = T)

vis_run(sta_day = "1982-01-01",
        end_day = "1984-12-31")

vis_run(sta_day = "2011-01-01",
        end_day = "2013-12-31")

dev.off()



#flux_stat----

#General

nc_flux_file <- paste0(run_dir, "output/mHM_Fluxes_States.nc")

nc_flux <- nc_open(nc_flux_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-02")
count_date <- length(date) - sta_date_ind


#Fluxes and states

snow_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
pet_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "PET")
aet_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "aET")
qto_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "Q")
sw1_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L01")
sw2_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L02")
sw3_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L03")
sw4_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L04")
sw5_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L05")
sw6_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SWC_L06")
pef_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "preEffect")
qdi_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "QD")
qba_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "QB")
qif_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "QIf")
qis_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "QIs")

snow_max <- apply(snow_cube, c(1, 2), max_na) #[mm]
pet_mea <- apply(pet_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
aet_mea <- apply(aet_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
qto_mea <- apply(qto_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
sw1_mea <- apply(sw1_cube, c(1, 2), mea_na) 
sw2_mea <- apply(sw2_cube, c(1, 2), mea_na) 
sw3_mea <- apply(sw3_cube, c(1, 2), mea_na) 
sw4_mea <- apply(sw4_cube, c(1, 2), mea_na) 
sw5_mea <- apply(sw5_cube, c(1, 2), mea_na) 
sw6_mea <- apply(sw6_cube, c(1, 2), mea_na) 
pef_mea <- apply(pef_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
qdi_mea <- apply(qdi_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
qba_mea <- apply(qba_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
qif_mea <- apply(qif_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
qis_mea <- apply(qis_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]

snow_max_c <- c(snow_max)
pet_mea_c <- c(pet_mea)
aet_mea_c <- c(aet_mea)
qto_mea_c <- c(qto_mea)
sw1_mea_c <- c(sw1_mea)
sw2_mea_c <- c(sw2_mea)
sw3_mea_c <- c(sw3_mea)
sw4_mea_c <- c(sw4_mea)
sw5_mea_c <- c(sw5_mea)
sw6_mea_c <- c(sw6_mea)
swc_mea_c <- sw1_mea_c + sw2_mea_c + sw3_mea_c + sw4_mea_c + sw5_mea_c + sw6_mea_c
pef_mea_c <- c(pef_mea)
qdi_mea_c <- c(qdi_mea)
qba_mea_c <- c(qba_mea)
qif_mea_c <- c(qif_mea)
qis_mea_c <- c(qis_mea)

range_et <- range(c(pet_mea_c, aet_mea_c), na.rm = T)

sn_tow_1 <- which(snow_max_c > 2500)
# sn_tow_2 <- which(snow_max_c > 2500)

snow_max_c[sn_tow_1] <- NA

cols_spat_sno <- foreach(i = 1:length(snow_max_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = snow_max_c[i],
          dat_ref = snow_max_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}

cols_spat_sno[sn_tow_1] <- "firebrick1"
# cols_spat_sno[sn_tow_2] <- "firebrick1"

cols_spat_pet <- foreach(i = 1:length(pet_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = pet_mea_c[i],
          dat_ref = range_et,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_aet <- foreach(i = 1:length(aet_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = aet_mea_c[i],
          dat_ref = range_et,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_qto <- foreach(i = 1:length(qto_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qto_mea_c[i],
          dat_ref = qto_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_swc <- foreach(i = 1:length(swc_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = swc_mea_c[i],
          dat_ref = swc_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_pef <- foreach(i = 1:length(pef_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = pef_mea_c[i],
          dat_ref = pef_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_qdi <- foreach(i = 1:length(qdi_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qdi_mea_c[i],
          dat_ref = qdi_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_qba <- foreach(i = 1:length(qba_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qba_mea_c[i],
          dat_ref = qba_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_qif <- foreach(i = 1:length(qif_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qif_mea_c[i],
          dat_ref = qif_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}
cols_spat_qis <- foreach(i = 1:length(qis_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qis_mea_c[i],
          dat_ref = qis_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}


pdf(paste0(bas_dir,"res_figs/flu_map.pdf"), width = 16, height = 18)

layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4,  rep(5, 7), 6,
                rep(7, 7), 8, rep(9, 7), 10, rep (11, 7), 12, 
                rep(13, 7), 14, rep(15, 7), 16, rep(17, 7), 18),
              3, 24, byrow = T), widths=c(), heights=c())
# layout.show(n=18)

par(family = "serif")

cex_pch <- 1.20
mar_1 <- c(1.5, 0.5, 1.5, 0.5)

#Plot: Snowpack
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_sno, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("a) Snow depth", side = 3, line = -1.0, cex = 1.5)

# st_sel_ind <- c(11502, 11505, 11508)
# points(c(lon)[st_sel_ind], c(lat)[st_sel_ind], pch = 4)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
# my_bre <- seq(range(log(snow_max_c), na.rm = T)[1], range(log(snow_max_c), na.rm = T)[2], length.out = length(my_col)+1)
my_bre <- seq(range(snow_max_c, na.rm = T)[1], range(snow_max_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(snow_max_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.50, 0), at = log(c(1, 10, 100, 1000, 2000)), labels = c(1, 10, 100, 1000, 2000), tck = -0.1, cex.axis = 1.6)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: PET
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_pet, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("b) pot. Evapotranspiration", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range_et[1], range_et[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(pet_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()


#Plot: aET
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_aet, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("c) act. Evapotranspiration", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range_et[1], range_et[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(pet_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()


#Plot: qto - total discharge generated per cell
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qto, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("d) Q total generated", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qto_mea_c, na.rm = T)[1], range(qto_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qto_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Effective precipitation
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_pef, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("e) effective Precipitation", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(pef_mea_c, na.rm = T)[1], range(pef_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(pef_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Direkt runoff
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qdi, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("f) Direct runoff", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qdi_mea_c, na.rm = T)[1], range(qdi_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qdi_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Baseflow generated
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qba, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("g) Baseflow generared", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qba_mea_c, na.rm = T)[1], range(qba_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qba_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Interflow fast
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qif, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("h) Interflow fast", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qif_mea_c, na.rm = T)[1], range(qif_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qif_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Interflow slow
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qis, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("i) Interflow slow", side = 3, line = -1.0, cex = 1.5)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qis_mea_c, na.rm = T)[1], range(qis_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qis_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

dev.off()


#Plot snow towers

st_sel_ind <- c(11502, 11505, 11508)

sel_x_1 <- NULL
sel_x_2 <- NULL
sel_x_3 <- NULL
for(i in 1:ncol(lon)){
  
  sel_dummy <- which(lon[, i] == c(lon)[st_sel_ind[1]])
  sel_x_1 <- c(sel_x_1, sel_dummy)
  
  sel_dummy <- which(lon[, i] == c(lon)[st_sel_ind[2]])
  sel_x_2 <- c(sel_x_2, sel_dummy)
  
  sel_dummy <- which(lon[, i] == c(lon)[st_sel_ind[3]])
  sel_x_3 <- c(sel_x_3, sel_dummy)
  
}

sel_y_1 <- NULL
sel_y_2 <- NULL
sel_y_3 <- NULL
for(i in 1:nrow(lon)){
  
  sel_dummy <- which(lon[i, ] == c(lon)[st_sel_ind[1]])
  sel_y_1 <- c(sel_y_1, sel_dummy)
  
  sel_dummy <- which(lon[i, ] == c(lon)[st_sel_ind[2]])
  sel_y_2 <- c(sel_y_2, sel_dummy)
  
  sel_dummy <- which(lon[i, ] == c(lon)[st_sel_ind[3]])
  sel_y_3 <- c(sel_y_3, sel_dummy)
  
}

stow_ts_1 <- ncvar_get(nc_flux, start = c(sel_x_1, sel_y_1, sta_date_ind), 
                       count = c(1, 1, count_date), varid = "snowpack")

stow_ts_2 <- ncvar_get(nc_flux, start = c(sel_x_2, sel_y_2, sta_date_ind), 
                       count = c(1, 1, count_date), varid = "snowpack")

stow_ts_3 <- ncvar_get(nc_flux, start = c(sel_x_3, sel_y_3, sta_date_ind), 
                       count = c(1, 1, count_date), varid = "snowpack")


pdf(paste0(bas_dir,"res_figs/snow_tower.pdf"), width = 12, height = 6)

par(family = "serif")
par(mfrow = c(3, 1))
par(mar = c(2, 4, 3, 1))

plot(date[-1], stow_ts_1, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[1]], 3), "  Lon: ", round(c(lon)[st_sel_ind[1]], 3)),
      side = 3, line = 0.2, cex = 1.2)

plot(date[-1], stow_ts_3, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[3]], 3), "  Lon: ", round(c(lon)[st_sel_ind[3]], 3)),
      side = 3, line = 0.2, cex = 1.2)

plot(date[-1], stow_ts_2, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[2]], 3), "  Lon: ", round(c(lon)[st_sel_ind[2]], 3)),
      side = 3, line = 0.2, cex = 1.2)

dev.off()

#seas_flod----

#get output fluxes from nc-file
nc_flux_file <- paste0(run_dir, "output/mHM_Fluxes_States.nc")

nc_flux <- nc_open(nc_flux_file)

lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-02")
count_date <- length(date) - sta_date_ind

snow_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
pef_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "preEffect")
qto_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "Q")


snow_max <- apply(snow_cube, c(1, 2), max_na) #[mm]
pef_mea <- apply(pef_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]
snow_max_c <- c(snow_max)
pef_mea_c <- c(pef_mea)

sn_tow_1 <- which(snow_max_c > 2500)
snow_max_c[sn_tow_1] <- NA

cols_spat_sno <- foreach(i = 1:length(snow_max_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = snow_max_c[i],
          dat_ref = snow_max_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}

cols_spat_sno[sn_tow_1] <- "firebrick1"

cols_spat_pef <- foreach(i = 1:length(pef_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = pef_mea_c[i],
          dat_ref = pef_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}

par(family = "serif")
cex_pch <- 1.20
mar_1 <- c(1.5, 0.5, 1.5, 0.5)

#Plot: Snowpack
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_sno, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("a) Snow depth", side = 3, line = -1.0, cex = 1.5)

rect(xleft = 7.1, xright = 10.3, ybottom = 46.2, ytop = 47.5)
rect(xleft = 5.3, xright = 12.0, ybottom = 48.1, ytop = 50.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
# my_bre <- seq(range(log(snow_max_c), na.rm = T)[1], range(log(snow_max_c), na.rm = T)[2], length.out = length(my_col)+1)
my_bre <- seq(range(snow_max_c, na.rm = T)[1], range(snow_max_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(snow_max_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.50, 0), at = log(c(1, 10, 100, 1000, 2000)), labels = c(1, 10, 100, 1000, 2000), tck = -0.1, cex.axis = 1.6)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Effective precipitation
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_pef, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("e) effective Precipitation", side = 3, line = -1.0, cex = 1.5)

rect(xleft = 7.1, xright = 10.3, ybottom = 46.2, ytop = 47.5)
rect(xleft = 5.3, xright = 12.0, ybottom = 48.1, ytop = 50.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(pef_mea_c, na.rm = T)[1], range(pef_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(pef_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Selcet cells
lat_in_1 <- c(lat)[which(c(lat) < 47.5 & c(lon) > 7.1)]
lat_in_2 <- c(lat)[which(c(lat) < 50.7 & c(lat) > 48.1)]

my_get_cube_col <- function(val_in, lats_in = lat, col_or_row = "col"){
  
  get_cube_index_col(val_in = val_in, lons_in = lats_in, col_or_row = col_or_row)
  
  }

my_get_cube_row <- function(val_in, lats_in = lat, col_or_row = "row"){
  
  get_cube_index_col(val_in = val_in, lons_in = lats_in, col_or_row = col_or_row)
  
}

#get index in cube from points inside sub-basin
cube_index_col_1 <- sapply(lat_in_1, my_get_cube_col)
cube_index_row_1 <- sapply(lat_in_1, my_get_cube_row)
cube_index_col_2 <- sapply(lat_in_2, my_get_cube_col)
cube_index_row_2 <- sapply(lat_in_2, my_get_cube_row)

#get effective precipitation + snow + runoff generated for nival part
for (i in 1:length(cube_index_col_1)) {
  
  print(paste(i, "of", length(cube_index_col_1)))
  
  epn_sing <- pef_cube [cube_index_col_1[i], cube_index_row_1[i], ]
  sno_sing <- snow_cube[cube_index_col_1[i], cube_index_row_1[i], ]
  qto_sing <- qto_cube [cube_index_col_1[i], cube_index_row_1[i], ]
  
  if(i == 1){
    epns <- epn_sing
    snos <- sno_sing
    qtos <- qto_sing
  }else{
    epns <- cbind(epns, epn_sing)
    snos <- cbind(snos, sno_sing)
    qtos <- cbind(qtos, qto_sing)
  }
  
}

niv_ep_sum <- apply(epns, 1, sum_na)
niv_qt_sum <- apply(qtos, 1, sum_na)
niv_sd_sum <- apply(snos, 1, sum_na)
niv_sd_sum_dif <- c(NA, diff(niv_sd_sum))
niv_sd_sum_dif[which(niv_sd_sum_dif > 0)] <- NA
niv_sn_sum <- niv_sd_sum_dif * -1 #melt positive values

#get effective precipitation + runoff generated for pluvial part
for (i in 1:length(cube_index_col_2)) {
  
  print(paste(i, "of", length(cube_index_col_2)))
  
  epn_sing <- pef_cube [cube_index_col_2[i], cube_index_row_2[i], ]
  qto_sing <- qto_cube [cube_index_col_2[i], cube_index_row_2[i], ]
  
  if(i == 1){
    epns_plu <- epn_sing
    qtos_plu <- qto_sing
  }else{
    epns_plu <- cbind(epns_plu, epn_sing)
    qtos_plu <- cbind(qtos_plu, qto_sing)
  }
}

plu_ep_sum <- apply(epns_plu, 1, sum_na)
plu_qt_sum <- apply(qtos_plu, 1, sum_na)


#Moving average sum
window_niv_sn <- 14
niv_ep_sum_ma <- rollapply(data = niv_ep_sum, width = window_niv_sn,
                           FUN = sum_na, align = "center", fill = NA)
window_niv_ep <- 14
niv_sn_sum_ma <- rollapply(data = niv_sn_sum, width = window_niv_ep,
                           FUN = sum_na, align = "center", fill = NA)
window_niv_qt <- 14
niv_qt_sum_ma <- rollapply(data = niv_qt_sum, width = window_niv_qt,
                           FUN = sum_na, align = "center", fill = NA)

plot(niv_ep_sum_ma, type = "l")
plot(niv_sn_sum_ma, type = "l")
plot(niv_qt_sum_ma, type = "l")

window_plu_ep <- 5
plu_ep_sum_ma <- rollapply(data = plu_ep_sum, width = window_plu_ep,
                        FUN = sum_na, align = "center", fill = NA)
window_plu_qt <- 5
plu_qt_sum_ma <- rollapply(data = plu_qt_sum, width = window_plu_qt,
                           FUN = sum_na, align = "center", fill = NA)

plot(plu_ep_sum_ma, type = "l")
plot(plu_qt_sum_ma, type = "l")

#Peak over threshold
pot_thre_niv_sn <- quantile(niv_sn_sum_ma, 0.95, na.rm = T)
pot_thre_niv_ep <- quantile(niv_ep_sum_ma, 0.95, na.rm = T)
pot_thre_niv_qt <- quantile(niv_qt_sum_ma, 0.95, na.rm = T)
pot_thre_plu_ep <- quantile(plu_ep_sum_ma, 0.95, na.rm = T)
pot_thre_plu_qt <- quantile(plu_qt_sum_ma, 0.95, na.rm = T)
# pot_thre_niv_ep <- pot_thre_niv_sn
# pot_thre_plu_ep <- pot_thre_niv_sn

pot_data_niv_sn <- data.frame(obs = niv_sn_sum_ma,
                              time = date[-1])
pot_data_niv_ep <- data.frame(obs = niv_ep_sum_ma,
                              time = date[-1])
pot_data_niv_qt <- data.frame(obs = niv_qt_sum_ma,
                              time = date[-1])
pot_data_plu_ep <- data.frame(obs = plu_ep_sum_ma,
                              time = date[-1])
pot_data_plu_qt <- data.frame(obs = plu_qt_sum_ma,
                              time = date[-1])

pot_data_niv_sn$obs[is.na(pot_data_niv_sn$obs)] <- 0
pot_data_niv_ep$obs[is.na(pot_data_niv_ep$obs)] <- 0
pot_data_niv_qt$obs[is.na(pot_data_niv_qt$obs)] <- 0
pot_data_plu_ep$obs[is.na(pot_data_plu_ep$obs)] <- 0
pot_data_plu_qt$obs[is.na(pot_data_plu_qt$obs)] <- 0

pot_peaks_niv_sn <- clust(data = pot_data_niv_sn, u = pot_thre_niv_sn, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_niv_ep <- clust(data = pot_data_niv_ep, u = pot_thre_niv_ep, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_niv_qt <- clust(data = pot_data_niv_qt, u = pot_thre_niv_qt, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_plu_ep <- clust(data = pot_data_plu_ep, u = pot_thre_plu_ep, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_plu_qt <- clust(data = pot_data_plu_qt, u = pot_thre_plu_qt, tim.cond = 14, clust.max = T, plot = F)

plot(niv_sn_sum_ma, type = "l")
points(pot_peaks_niv_sn[, 3], pot_peaks_niv_sn[, 2], col = "red3")

plot(niv_ep_sum_ma, type = "l")
points(pot_peaks_niv_ep[, 3], pot_peaks_niv_ep[, 2], col = "orange3")

plot(niv_qt_sum_ma, type = "l")
points(pot_peaks_niv_qt[, 3], pot_peaks_niv_qt[, 2], col = "gold3")

plot(plu_ep_sum_ma, type = "l")
points(pot_peaks_plu_ep[, 3], pot_peaks_plu_ep[, 2], col = "blue3")

plot(plu_qt_sum_ma, type = "l")
points(pot_peaks_plu_qt[, 3], pot_peaks_plu_qt[, 2], col = "steelblue3")

peaks_doy_niv_sn <- as.numeric(format(pot_data_niv_sn$time[pot_peaks_niv_sn[, 3]], '%j'))
peaks_doy_niv_ep <- as.numeric(format(pot_data_niv_ep$time[pot_peaks_niv_ep[, 3]], '%j'))
peaks_doy_niv_qt <- as.numeric(format(pot_data_niv_qt$time[pot_peaks_niv_qt[, 3]], '%j'))
peaks_doy_plu_ep <- as.numeric(format(pot_data_plu_ep$time[pot_peaks_plu_ep[, 3]], '%j'))
peaks_doy_plu_qt <- as.numeric(format(pot_data_plu_qt$time[pot_peaks_plu_qt[, 3]], '%j'))

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

ylims <- range(c(pot_peaks_niv_sn[ ,2], pot_peaks_niv_ep[ ,2], pot_peaks_niv_qt[ ,2], 
                 pot_peaks_plu_ep[, 2], pot_peaks_plu_qt[, 2]))
plot(peaks_doy_niv_sn, pot_peaks_niv_sn[, 2], xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims, type = "n")
points(peaks_doy_niv_sn, pot_peaks_niv_sn[, 2], col = alpha("darkred", alpha = 0.3), pch = 19, cex = 2)
points(peaks_doy_niv_ep, pot_peaks_niv_ep[, 2], col = alpha("orange3", alpha = 0.3), pch = 19, cex = 2)
points(peaks_doy_niv_qt, pot_peaks_niv_qt[, 2], col = alpha("gold3", alpha = 0.3), pch = 19, cex = 2)
points(peaks_doy_plu_ep, pot_peaks_plu_ep[, 2], col = alpha("darkblue", alpha = 0.3), pch = 19, cex = 2)
points(peaks_doy_plu_qt, pot_peaks_plu_qt[, 2], col = alpha("steelblue4", alpha = 0.3), pch = 19, cex = 2)

axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
box()

d_1 <- density(peaks_doy_niv_sn)
d_2 <- density(peaks_doy_plu_qt)

plot(d_1, xlim = c(0, 1000), ylim = c(0,0.01))
lines(d_2)
d <- density(peaks_doy_plu_qt)
plot(d, xlim = c(0, 365))

# #Annual maxima
# data_day_niv_ep <- ord_day(data_in = niv_ep_sum_ma,
#                            date = date[-1],
#                            start_y = 1954,
#                            end_y = 2013,
#                            break_day = 0,
#                            do_ma = F,
#                            window_width = 30)
# 
# data_day_niv_sn <- ord_day(data_in = niv_sn_sum_ma,
#                            date = date[-1],
#                            start_y = 1954,
#                            end_y = 2013,
#                            break_day = 0,
#                            do_ma = F,
#                            window_width = 30)
# 
# data_day_plu_ep <- ord_day(data_in = plu_ep_sum_ma,
#                            date = date[-1],
#                            start_y = 1954,
#                            end_y = 2013,
#                            break_day = 0,
#                            do_ma = F,
#                            window_width = 30)
# 
# yea_mag_niv_ep <- apply(data_day_niv_ep, 1, max_na)
# yea_mag_niv_sn <- apply(data_day_niv_sn, 1, max_na)
# yea_mag_plu_ep <- apply(data_day_plu_ep, 1, max_na)
# 
# max_doy <- function(data_in){
#   
#   doy_max <- which(data_in == max_na(data_in))[1]
#   
#   return(doy_max)
# }
# 
# yea_doy_niv_ep <- apply(data_day_niv_ep, 1, max_doy)
# yea_doy_niv_sn <- apply(data_day_niv_sn, 1, max_doy)
# yea_doy_plu_ep <- apply(data_day_plu_ep, 1, max_doy)
# 
# plot(yea_doy_niv_sn, yea_mag_niv_sn, xlim = c(0, 365), col = "red3", pch = 19)
# 
# plot(yea_doy_niv_ep, yea_mag_niv_ep, xlim = c(0, 365), col = "orange3", pch = 19)
# 
# plot(yea_doy_plu_ep, yea_mag_plu_ep, xlim = c(0, 365), col = "blue3", pch = 19)