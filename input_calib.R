###

#Analyze mhm model results

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales, ncdf4,
               ncdf4.helpers, sp, raster, viridis, meltimr, POT, readr, hydroGOF,
               CoinCalc, seas)

# run_dir <- "D:/nrc_user/rottler/mhm_run/6935053/"
run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"
out_dir <- "D:/nrc_user/rottler/mhm_run/6435060/output/OBS/"

bas_dir <- "U:/rhine_fut/R/"

#directory to snow cover data from DLR
scf_dlr_dir <- "D:/nrc_user/rottler/SCF_data/snow_dlr/SnowPack_DLR.tar/SnowPack_DLR/" 

#load functions
source(paste0(bas_dir, "mhm_rhine/functs.R"))

sta_yea <- 1951
end_yea <- 2013

date_simu <- seq(as.Date("1951-01-01", format = "%Y-%m-%d"), 
                 as.Date("2013-12-31", format = "%Y-%m-%d"), by = "day")

stopCluster(my_clust)

n_cores <- 25 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#inp_vis----

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
  
  val2col(val_in = c(temps_mea)[i],
          dat_ref = c(temps_mea),
          do_bicol = F,
          virid_dir = 1)
  
}
cols_spat_pre <- foreach(i = 1:length(c(precs_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(precs_mea)[i],
          dat_ref = c(precs_mea),
          do_bicol = F)
  
}
cols_spat_eva <- foreach(i = 1:length(c(evapo_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(evapo_mea)[i],
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

plot(facc_bas, col = c("grey92", "blue2", "darkblue"), breaks = c(0, 5000, 1000000, 10000000), axes=F, legend = T, ylab = "", xlab = "", box = F)
mtext("Flow accumulation", side = 3, line = 0.0)

dev.off()



#dis_vis----

dis_mhm <- read.table(paste0(out_dir, "daily_discharge.out"), header = T)

NSE(dis_mhm$Qsim_0006435060, dis_mhm$Qobs_0006435060)
NSE(dis_mhm$Qsim_0006935051, dis_mhm$Qobs_0006935051)


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
# sn_tow_1
st_sel_ind <- c(11502, 11319, 11506)

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

plot(date[-1], stow_ts_2, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[2]], 3), "  Lon: ", round(c(lon)[st_sel_ind[2]], 3)),
      side = 3, line = 0.2, cex = 1.2)

plot(date[-1], stow_ts_3, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[3]], 3), "  Lon: ", round(c(lon)[st_sel_ind[3]], 3)),
      side = 3, line = 0.2, cex = 1.2)

dev.off()

#dis_rout----

nc_disc_file <- paste0(out_dir, "mRM_Fluxes_States.nc")
nc_disc <- nc_open(nc_disc_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_disc, varid = "lon")
lat <- ncdf4::ncvar_get(nc_disc, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_disc, time.dim.name = "time")))

count_date <- length(date)

disc_cube <- ncvar_get(nc_disc, start = c(1, 1, 1), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "Qrouted")

disc_mea <- apply(disc_cube, c(1,2), mea_na)

#Get meta data and measured time series for selected gauges (GRDC)

grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

lobi_file <- paste0(grdc_dir, "6435060_Q_Day.Cmd.txt")
koel_file <- paste0(grdc_dir, "6335060_Q_Day.Cmd.txt")
coch_file <- paste0(grdc_dir, "6336050_Q_Day.Cmd.txt")
kaub_file <- paste0(grdc_dir, "6335100_Q_Day.Cmd.txt")
wuer_file <- paste0(grdc_dir, "6335500_Q_Day.Cmd.txt")
worm_file <- paste0(grdc_dir, "6335180_Q_Day.Cmd.txt")
rock_file <- paste0(grdc_dir, "6335600_Q_Day.Cmd.txt")
spey_file <- paste0(grdc_dir, "6335170_Q_Day.Cmd.txt")
base_file <- paste0(grdc_dir, "6935051_Q_Day.Cmd.txt")
unte_file <- paste0(grdc_dir, "6935300_Q_Day.Cmd.txt")
reki_file <- paste0(grdc_dir, "6935054_Q_Day.Cmd.txt")

file_paths <- c(lobi_file, koel_file, coch_file, kaub_file, wuer_file, worm_file,
                rock_file, spey_file, base_file, unte_file, reki_file)

grdc_meta <- NULL

for(i in 1:length(file_paths)){
  
  #get rows with meta information
  meta_rows <- read_lines(file_paths[i], n_max = 32)
  meta_rows <- iconv(meta_rows, "UTF-8", "ASCII", "")
  
  #Name
  sta_name <- substr(meta_rows[11], 26, nchar(meta_rows[11]))
  
  #Longitude
  sta_long <- substr(meta_rows[14], 24, nchar( meta_rows[14]))
  
  #Latitude
  sta_lati <- substr(meta_rows[13], 24, nchar(meta_rows[13]))
  
  #Meta data single station
  meta_sing <- c(sta_name, sta_lati, sta_long)
  
  #Collect meta data all stations
  grdc_meta <- rbind(grdc_meta, meta_sing)
  
}

colnames(grdc_meta) <- c("name", "latitude", "longitude")
rownames(grdc_meta) <- NULL
grdc_meta <- as.data.frame(grdc_meta)
grdc_meta$name  <- as.character(levels(grdc_meta$name))[grdc_meta$name]
grdc_meta$latitude   <- as.numeric(levels(grdc_meta$latitude))[grdc_meta$latitude]
grdc_meta$longitude  <- as.numeric(levels(grdc_meta$longitude))[grdc_meta$longitude]

disc_lobi_full <- read_grdc(lobi_file)
disc_koel_full <- read_grdc(koel_file)
disc_coch_full <- read_grdc(coch_file)
disc_kaub_full <- read_grdc(kaub_file)
disc_wuer_full <- read_grdc(wuer_file)
disc_worm_full <- read_grdc(worm_file)
disc_rock_full <- read_grdc(rock_file)
disc_spey_full <- read_grdc(spey_file)
disc_base_full <- read_grdc(base_file)
disc_unte_full <- read_grdc(unte_file)
disc_reki_full <- read_grdc(reki_file)

disc_lobi <- disc_lobi_full[which(disc_lobi_full$date %in% date_simu), ]
disc_koel <- disc_koel_full[which(disc_koel_full$date %in% date_simu), ]
disc_coch <- disc_coch_full[which(disc_coch_full$date %in% date_simu), ]
disc_kaub <- disc_kaub_full[which(disc_kaub_full$date %in% date_simu), ]
disc_wuer <- disc_wuer_full[which(disc_wuer_full$date %in% date_simu), ]
disc_worm <- disc_worm_full[which(disc_worm_full$date %in% date_simu), ]
disc_rock <- disc_rock_full[which(disc_rock_full$date %in% date_simu), ]
disc_spey <- disc_spey_full[which(disc_spey_full$date %in% date_simu), ]
disc_base <- disc_base_full[which(disc_base_full$date %in% date_simu), ]
disc_unte <- disc_unte_full[which(disc_unte_full$date %in% date_simu), ]
disc_reki <- disc_reki_full[which(disc_reki_full$date %in% date_simu), ]

#Get simulated runoff for selected gauges

gaug_spa <- SpatialPoints(grdc_meta[, c(2, 3)], proj4string = crswgs84)

grid_spa <- SpatialPoints(cbind(c(lat), c(lon)), proj4string = crswgs84)

#determine cells (rows-clumns in grid) representing gauges
f_index_row <- function (val_in, lons_in = lon, col_or_row = "row"){
  if (col_or_row == "row") {
    index_out <- which(round(lons_in, digits = 6) == round(val_in, 
                                                           digits = 6), arr.ind = T)[1, 1]
  }
  if (col_or_row == "col") {
    index_out <- which(round(lons_in, digits = 6) == round(val_in, 
                                                           digits = 6), arr.ind = T)[1, 2]
  }
  return(index_out)
}
f_index_col <- function (val_in, lons_in = lon, col_or_row = "col"){
  if (col_or_row == "row") {
    index_out <- which(round(lons_in, digits = 6) == round(val_in, 
                                                           digits = 6), arr.ind = T)[1, 1]
  }
  if (col_or_row == "col") {
    index_out <- which(round(lons_in, digits = 6) == round(val_in, 
                                                           digits = 6), arr.ind = T)[1, 2]
  }
  return(index_out)
}

coords_sel_gaugs <- NULL
rows_sel_gaugs <- NULL
cols_sel_gaugs <- NULL

for(i in 1:length(gaug_spa)){

  cell_sel <- which(pointDistance(gaug_spa@coords[i, c(2, 1)], grid_spa@coords[, c(2, 1)], lonlat = T) ==
                      min_na(pointDistance(gaug_spa@coords[i, c(2, 1)], grid_spa@coords[, c(2, 1)], lonlat = T)))
  
  row_sel <- f_index_row(grid_spa@coords[cell_sel, 2])
  col_sel <- f_index_col(grid_spa@coords[cell_sel, 2])
  
  coords_sel_gaugs <- rbind(coords_sel_gaugs, grid_spa@coords[cell_sel, ])
  rows_sel_gaugs <- c(rows_sel_gaugs, row_sel)
  cols_sel_gaugs <- c(cols_sel_gaugs, col_sel)
  
}

#gauge cochem one row lower
rows_sel_gaugs[3] <- rows_sel_gaugs[3]+1

simu_lobi <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[1], cols_sel_gaugs[1], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_koel <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[2], cols_sel_gaugs[2], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_coch <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[3], cols_sel_gaugs[3], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_kaub <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[4], cols_sel_gaugs[4], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_wuer <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[5], cols_sel_gaugs[5], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_worm <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[6], cols_sel_gaugs[6], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_rock <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[7], cols_sel_gaugs[7], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_spey <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[8], cols_sel_gaugs[8], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_base <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[9], cols_sel_gaugs[9], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_unte <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[10], cols_sel_gaugs[10], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")
simu_reki <- ncvar_get(nc_disc, start = c(rows_sel_gaugs[11], cols_sel_gaugs[11], 1), 
                       count = c(1, 1, count_date), varid = "Qrouted")

#Validation different time windows
sta_yea_val_all <- c(1954, 1954, 1974, 1994)
end_yea_val_all <- c(2013, 1973, 1993, 2013)

kge_out <- NULL
nse_out <- NULL
for(i in 1:length(sta_yea_val_all)){

  date_vali <- seq(as.Date(paste0(sta_yea_val_all[i],"-01", "-01")), 
                   as.Date(paste0(end_yea_val_all[i],"-12", "-31")), 
                   by = "day")
  date_sel_ind <- which(date_simu %in% date_vali)
  
  kge_sel <-  c(
    KGE(as.numeric(simu_lobi[date_sel_ind]), as.numeric(disc_lobi$value[date_sel_ind])),
    KGE(as.numeric(simu_koel[date_sel_ind]), as.numeric(disc_koel$value[date_sel_ind])),
    KGE(as.numeric(simu_coch[date_sel_ind]), as.numeric(disc_coch$value[date_sel_ind])),
    KGE(as.numeric(simu_kaub[date_sel_ind]), as.numeric(disc_kaub$value[date_sel_ind])),
    KGE(as.numeric(simu_wuer[date_sel_ind]), as.numeric(disc_wuer$value[date_sel_ind])),
    KGE(as.numeric(simu_worm[date_sel_ind]), as.numeric(disc_worm$value[date_sel_ind])),
    KGE(as.numeric(simu_rock[date_sel_ind]), as.numeric(disc_rock$value[date_sel_ind])),
    KGE(as.numeric(simu_spey[date_sel_ind]), as.numeric(disc_spey$value[date_sel_ind])),
    KGE(as.numeric(simu_base[date_sel_ind]), as.numeric(disc_base$value[date_sel_ind])),
    KGE(as.numeric(simu_unte[date_sel_ind]), as.numeric(disc_unte$value[date_sel_ind])),
    KGE(as.numeric(simu_reki[date_sel_ind]), as.numeric(disc_reki$value[date_sel_ind]))
  )  

  nse_sel <-  c(
    NSE(as.numeric(simu_lobi[date_sel_ind]), as.numeric(disc_lobi$value[date_sel_ind])),
    NSE(as.numeric(simu_koel[date_sel_ind]), as.numeric(disc_koel$value[date_sel_ind])),
    NSE(as.numeric(simu_coch[date_sel_ind]), as.numeric(disc_coch$value[date_sel_ind])),
    NSE(as.numeric(simu_kaub[date_sel_ind]), as.numeric(disc_kaub$value[date_sel_ind])),
    NSE(as.numeric(simu_wuer[date_sel_ind]), as.numeric(disc_wuer$value[date_sel_ind])),
    NSE(as.numeric(simu_worm[date_sel_ind]), as.numeric(disc_worm$value[date_sel_ind])),
    NSE(as.numeric(simu_rock[date_sel_ind]), as.numeric(disc_rock$value[date_sel_ind])),
    NSE(as.numeric(simu_spey[date_sel_ind]), as.numeric(disc_spey$value[date_sel_ind])),
    NSE(as.numeric(simu_base[date_sel_ind]), as.numeric(disc_base$value[date_sel_ind])),
    NSE(as.numeric(simu_unte[date_sel_ind]), as.numeric(disc_unte$value[date_sel_ind])),
    NSE(as.numeric(simu_reki[date_sel_ind]), as.numeric(disc_reki$value[date_sel_ind]))
  ) 
  
  kge_out <- cbind(kge_out, kge_sel)
  nse_out <- cbind(nse_out, nse_sel)
}

round(kge_out, 3)
round(nse_out, 3)

#Read basins and river network
basin_lobi_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/lobith_catch.shp")
basin_koel_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/koeln_catch.shp")
basin_coch_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/cochem_catch.shp")
basin_kaub_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/kaub_catch.shp")
basin_wuer_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/wuerzburg_catch.shp")
basin_worm_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/worms_catch.shp")
basin_rock_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/rockenau_catch.shp")
basin_spey_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/speyer_catch.shp")
basin_base_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/basel_catch.shp")
basin_unte_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/untersiggenthal_catch.shp")
basin_reki_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/rekingen_catch.shp")
river_netw_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/river_network.shp")

basin_lobi <- spTransform(basin_lobi_raw, CRS = crswgs84)
basin_koel <- spTransform(basin_koel_raw, CRS = crswgs84)
basin_coch <- spTransform(basin_coch_raw, CRS = crswgs84)
basin_kaub <- spTransform(basin_kaub_raw, CRS = crswgs84)
basin_wuer <- spTransform(basin_wuer_raw, CRS = crswgs84)
basin_worm <- spTransform(basin_worm_raw, CRS = crswgs84)
basin_rock <- spTransform(basin_rock_raw, CRS = crswgs84)
basin_spey <- spTransform(basin_spey_raw, CRS = crswgs84)
basin_base <- spTransform(basin_base_raw, CRS = crswgs84)
basin_unte <- spTransform(basin_unte_raw, CRS = crswgs84)
basin_reki <- spTransform(basin_reki_raw, CRS = crswgs84)
river_netw <- spTransform(river_netw_raw, CRS = crswgs84)

cols_spat_dis <- foreach(i = 1:length(c(disc_mea)), .combine = 'cbind') %dopar% {
  
  val2col(val_in = c(disc_mea)[i],
          dat_ref = c(disc_mea),
          do_log = T,
          do_bicol = F,
          virid_dir = -1)
  
}


pdf(paste0(bas_dir,"res_figs/dis_vali_raw.pdf"), width = 12, height = 6)

layout(matrix(c(rep(1, 7), 2, rep(3, 8)),
              1, 16, byrow = T), widths=c(), heights=c())

#Plot routed discharge

par(family = "serif")

cex_pch <- 1.32

par(mar = c(0.5, 0.5, 1.0, 0.5))

plot(c(lon), c(lat), pch = 15, col = cols_spat_dis, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("a) Discharge routed", side = 3, line = -1.0, cex = 1.5)
points(coords_sel_gaugs[, 2], coords_sel_gaugs[, 1], pch = 25, bg = "white", cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))

my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
# my_bre <- seq(range(log(snow_max_c), na.rm = T)[1], range(log(snow_max_c), na.rm = T)[2], length.out = length(my_col)+1)
my_bre <- seq(range(log(c(disc_mea)), na.rm = T)[1], range(log(c(disc_mea)), na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(c(disc_mea)), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.50, 0), at = log(c(1, 10, 100, 1000, 2000)), labels = c(1, 10, 100, 1000, 2000), tck = -0.1, cex.axis = 1.6)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[ln(m"^"3", "s"^"-1",")]")), side = 3, line = 0.7, cex = 1.2)
box()

#Plot Basins + River network

par(mar = c(0.5, 0.5, 1.0, 0.5))

col_rhine <- alpha("black", alpha = 0.15)
col_tribu <- alpha("steelblue4", alpha = 0.5)

plot(basin_lobi, col = col_rhine, border = F)
plot(basin_koel, col = col_rhine, border = F, add = T)
plot(basin_kaub, col = col_rhine, border = F, add = T)
plot(basin_worm, col = col_rhine, border = F, add = T)
plot(basin_spey, col = col_rhine, border = F, add = T)
plot(basin_base, col = col_rhine, border = F, add = T)
plot(basin_reki, col = col_rhine, border = F, add = T)

plot(basin_coch, col = col_tribu, border = F, add = T)
plot(basin_wuer, col = col_tribu, border = F, add = T)
plot(basin_rock, col = col_tribu, border = F, add = T)
plot(basin_unte, col = col_tribu, border = F, add = T)

plot(river_netw, col = "darkblue", add = T)

points(coords_sel_gaugs[, 2], coords_sel_gaugs[, 1], pch = 25, bg = "black", cex = 1.7)
points(coords_sel_gaugs[c(3, 5, 7, 10), 2], coords_sel_gaugs[c(3, 5, 7, 10), 1], pch = 25, col = "steelblue4", bg = "steelblue4", cex = 1.7)
mtext("b) Gauges validation", side = 3, line = -1.0, cex = 1.5)

dev.off()


#disc_quan----

#Select gauge/time series

quants <- seq(0.01, 0.99, by = 0.01)
date_simu <- seq(as.Date("1954-01-01", format = "%Y-%m-%d"), 
                 as.Date("2013-12-31", format = "%Y-%m-%d"), by = "day")

f_qvalu_obs <- function(quant_sel){dis_ana(disc = disc_obs_sel$value,
                                           date = disc_obs_sel$date,
                                           start_year = sta_yea,
                                           end_year = end_yea,
                                           break_day = 0,
                                           quant_in = quant_sel,
                                           do_moving_average = F,
                                           window_width = 30,
                                           method_analys = "quantile",
                                           method_quant = "empirical"
)}
f_qvalu_sim <- function(quant_sel){dis_ana(disc = disc_sim_sel,
                                           date = date_simu,
                                           start_year = sta_yea,
                                           end_year = end_yea,
                                           break_day = 0,
                                           quant_in = quant_sel,
                                           do_moving_average = F,
                                           window_width = 30,
                                           method_analys = "quantile",
                                           method_quant = "empirical"
)}

disc_obs_sel <- disc_reki
disc_sim_sel <- simu_reki

qvalu_obs_reki <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

qvalu_sim_reki <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif_reki <- qvalu_obs_reki - qvalu_sim_reki

disc_obs_sel <- disc_base
disc_sim_sel <- simu_base

qvalu_obs_base <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

qvalu_sim_base <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif_base <- qvalu_obs_base - qvalu_sim_base

disc_obs_sel <- disc_coch
disc_sim_sel <- simu_coch

qvalu_obs_coch <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

qvalu_sim_coch <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif_coch <- qvalu_obs_coch - qvalu_sim_coch

disc_obs_sel <- disc_rock
disc_sim_sel <- simu_rock

qvalu_obs_rock <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

qvalu_sim_rock <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif_rock <- qvalu_obs_rock - qvalu_sim_rock

disc_obs_sel <- disc_koel
disc_sim_sel <- simu_koel

qvalu_obs_koel <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_obs(k)
}

qvalu_sim_koel <- foreach(k = quants, .combine = 'cbind') %dopar%{
  f_qvalu_sim(k)
}

qvalu_dif_koel <- qvalu_obs_koel - qvalu_sim_koel


#Plot: Seasonality of runoff

pdf(paste0(bas_dir, "res_figs/runoff_quant",".pdf"), width = 2*6.0, height = 6*1.7)
# tiff(paste0(bas_dir, "res_figs/runoff_qu_", stat_sel,".tiff"), width = 3*2.0, height = 3*1.5,
#      units = "in", res = 800)
par(family = "serif")
 
layout(matrix(c(rep(1,  8),  2, rep(7,  8),  8,
                rep(3,  8),  4, rep(9,  8), 10,
                rep(5,  8),  6, rep(11, 8), 12,
                rep(25, 18),
                rep(13, 8), 14, rep(19, 8), 20,
                rep(15, 8), 16, rep(21, 8), 22,
                rep(17, 8), 18, rep(23, 8), 24),
              7, 18, byrow = T), widths=c(), heights=c(1, 1, 1, 0.2, 1, 1, 1))

#Basel Observations
stat_sel <- "Basel"
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_obs_base, qvalu_sim_base)), 
               alptempr::max_na(c(qvalu_obs_base, qvalu_sim_base)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_obs_base, cols = my_col, breaks = my_bre, 
          header = paste("a)",  stat_sel, "observations"), lab_unit = "[m³/s]")

#Basel Simulations
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_sim_base, qvalu_obs_base)), 
               alptempr::max_na(c(qvalu_sim_base, qvalu_obs_base)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_sim_base, cols = my_col, breaks = my_bre, 
          header = paste("b)", stat_sel, "simulations"), lab_unit = "[m³/s]")

#Basel Difference
cols_max <- colorRampPalette(c(rep("grey98", 15), "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,2,3,4)], "lightskyblue3", rep("lightcyan3", 1), rep("grey98", 15)))(100)
my_col <- c(cols_min, cols_max)
# my_bre <- seq(-max_na(abs(qvalu_dif)), max_na(abs(qvalu_dif)), length.out = length(my_col)+1)
my_bre <- c(-(lseq(0.01, max_na(abs(qvalu_dif_base)), length.out = length(my_col)/2)[(length(my_col)/2):1]),
            lseq(0.01, max_na(abs(qvalu_dif_base)), length.out = length(my_col)/2+1))

dis_image(data_plot = qvalu_dif_base, cols = my_col, breaks = my_bre, 
          header = paste("c)", stat_sel, "obs - sim"), lab_unit = "[m³/s]", do_cont = F)

#Rockenau Observations
stat_sel <- "Rockenau"
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_obs_rock, qvalu_sim_rock)), 
               alptempr::max_na(c(qvalu_obs_rock, qvalu_sim_rock)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_obs_rock, cols = my_col, breaks = my_bre, 
          header = paste("d)",  stat_sel, "observations"), lab_unit = "[m³/s]")

#Rockenau Simulations
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_sim_rock, qvalu_obs_rock)), 
               alptempr::max_na(c(qvalu_sim_rock, qvalu_obs_rock)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_sim_rock, cols = my_col, breaks = my_bre, 
          header = paste("e)", stat_sel, "simulations"), lab_unit = "[m³/s]")

#Rockenau Difference
cols_max <- colorRampPalette(c(rep("grey98", 15), "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,2,3,4)], "lightskyblue3", rep("lightcyan3", 1), rep("grey98", 15)))(100)
my_col <- c(cols_min, cols_max)
# my_bre <- seq(-max_na(abs(qvalu_dif)), max_na(abs(qvalu_dif)), length.out = length(my_col)+1)
my_bre <- c(-(lseq(0.01, max_na(abs(qvalu_dif_rock)), length.out = length(my_col)/2)[(length(my_col)/2):1]),
            lseq(0.01, max_na(abs(qvalu_dif_rock)), length.out = length(my_col)/2+1))

dis_image(data_plot = qvalu_dif_rock, cols = my_col, breaks = my_bre, 
          header = paste("f)", stat_sel, "obs - sim"), lab_unit = "[m³/s]", do_cont = F)

#Cochem Observations
stat_sel <- "Cochem"
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_obs_coch, qvalu_sim_coch)), 
               alptempr::max_na(c(qvalu_obs_coch, qvalu_sim_coch)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_obs_coch, cols = my_col, breaks = my_bre, 
          header = paste("g)",  stat_sel, "observations"), lab_unit = "[m³/s]")

#Cochem Simulations
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_sim_coch, qvalu_obs_coch)), 
               alptempr::max_na(c(qvalu_sim_coch, qvalu_obs_coch)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_sim_coch, cols = my_col, breaks = my_bre, 
          header = paste("h)", stat_sel, "simulations"), lab_unit = "[m³/s]")

#Cochem Difference
cols_max <- colorRampPalette(c(rep("grey98", 15), "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,2,3,4)], "lightskyblue3", rep("lightcyan3", 1), rep("grey98", 15)))(100)
my_col <- c(cols_min, cols_max)
# my_bre <- seq(-max_na(abs(qvalu_dif)), max_na(abs(qvalu_dif)), length.out = length(my_col)+1)
my_bre <- c(-(lseq(0.01, max_na(abs(qvalu_dif_coch)), length.out = length(my_col)/2)[(length(my_col)/2):1]),
            lseq(0.01, max_na(abs(qvalu_dif_coch)), length.out = length(my_col)/2+1))

dis_image(data_plot = qvalu_dif_coch, cols = my_col, breaks = my_bre, 
          header = paste("i)", stat_sel, "obs - sim"), lab_unit = "[m³/s]", do_cont = F)


#Cologne Observations
stat_sel <- "Cologne"
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_obs_koel, qvalu_sim_koel)), 
               alptempr::max_na(c(qvalu_obs_koel, qvalu_sim_koel)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_obs_koel, cols = my_col, breaks = my_bre, 
          header = paste("j)",  stat_sel, "observations"), lab_unit = "[m³/s]")

#Cologne Simulations
cols_max <- grDevices::colorRampPalette(c("white", "cadetblue3", viridis::viridis(9, direction = 1)[c(4:1, 1)]))(100)
cols_min <- grDevices::colorRampPalette(c("red4","orangered4", "orange2","gold2", "yellow2", "white"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- lseq(alptempr::min_na(c(qvalu_sim_koel, qvalu_obs_koel)), 
               alptempr::max_na(c(qvalu_sim_koel, qvalu_obs_koel)), length.out = length(my_col)+1)

dis_image(data_plot = qvalu_sim_koel, cols = my_col, breaks = my_bre, 
          header = paste("k)", stat_sel, "simulations"), lab_unit = "[m³/s]")

#Cologne Difference
cols_max <- colorRampPalette(c(rep("grey98", 15), "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,2,3,4)], "lightskyblue3", rep("lightcyan3", 1), rep("grey98", 15)))(100)
my_col <- c(cols_min, cols_max)
# my_bre <- seq(-max_na(abs(qvalu_dif)), max_na(abs(qvalu_dif)), length.out = length(my_col)+1)
my_bre <- c(-(lseq(0.01, max_na(abs(qvalu_dif_koel)), length.out = length(my_col)/2)[(length(my_col)/2):1]),
            lseq(0.01, max_na(abs(qvalu_dif_koel)), length.out = length(my_col)/2+1))

dis_image(data_plot = qvalu_dif_koel, cols = my_col, breaks = my_bre, 
          header = paste("l)", stat_sel, "obs - sim"), lab_unit = "[m³/s]", do_cont = F)

dev.off()


#seas_flood----

#Get basin .shp
basin_coch_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/cochem_catch.shp")
basin_base_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/basel_catch.shp")
basin_coch <- spTransform(basin_coch_raw, CRS = crswgs84)
basin_base <- spTransform(basin_base_raw, CRS = crswgs84)

#get simulation results from output nc-file
nc_flux_file <- paste0(run_dir, "output/mHM_Fluxes_States.nc")

nc_flux <- nc_open(nc_flux_file)

lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-02")
count_date <- length(date)

#snowpack
snow_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
#effective precipitation
pef_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "preEffect")
#total runoff generated
qto_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "Q")
#total runoff generated
qto_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "Q")

sw1_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L01")
sw2_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L02")
sw3_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L03")
sw4_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L04")
sw5_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L05")
sw6_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                      count = c(nrow(lon), ncol(lon), count_date), varid = "SM_L06")

#get precipitation input
nc_prec_file <- paste0(run_dir, "input/meteo/pre.nc")

nc_prec <- nc_open(nc_prec_file)

date <- as.Date(as.character(nc.get.time.series(nc_prec, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-01")
count_date <- length(date) - sta_date_ind

pre_cube <- ncvar_get(nc_prec, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lat), count_date), varid = "pre")

pre_cube <- pre_cube[ , , which(date_simu %in% date)] #selecte simulatin period

#spatial grid points from lat/lon
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)
# grid_points_cube     <- sp::spTransform(grid_points_cube_84, CRS = crs(epsg3035, asText = T))

#grid points inside watersheds
inside_base <- !is.na(sp::over(grid_points_cube_84, as(basin_base, "SpatialPolygons")))
inside_coch <- !is.na(sp::over(grid_points_cube_84, as(basin_coch, "SpatialPolygons")))
grid_points_base <- grid_points_cube_84[which(inside_base == T)]
grid_points_coch <- grid_points_cube_84[which(inside_coch == T)]

#Select cells for Basel/Cochem watershed
lat_in_base <- grid_points_base@coords[, 2]
lat_in_coch <- grid_points_coch@coords[, 2]

my_get_cube_col <- function(val_in, lats_in = lat, col_or_row = "col"){
  
  get_cube_index_col(val_in = val_in, lons_in = lats_in, col_or_row = col_or_row)
  
}

my_get_cube_row <- function(val_in, lats_in = lat, col_or_row = "row"){
  
  get_cube_index_col(val_in = val_in, lons_in = lats_in, col_or_row = col_or_row)
  
}

#get index in cube from points inside sub-basins
cube_index_col_base <- sapply(lat_in_base, my_get_cube_col)
cube_index_row_base <- sapply(lat_in_base, my_get_cube_row)
cube_index_col_coch <- sapply(lat_in_coch, my_get_cube_col)
cube_index_row_coch <- sapply(lat_in_coch, my_get_cube_row)

#get simulation results for Basel watershed
for (i in 1:length(cube_index_col_base)) {
  
  print(paste(i, "of", length(cube_index_col_base)))
  
  epn_sing <- pef_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  sno_sing <- snow_cube[cube_index_col_base[i], cube_index_row_base[i], ]
  qto_sing <- qto_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  pre_sing <- pre_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  sw1_sing <- sw1_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  sw2_sing <- sw2_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  sw3_sing <- sw3_cube [cube_index_col_base[i], cube_index_row_base[i], ]
  
  if(i == 1){
    epn_base <- epn_sing
    sno_base <- sno_sing
    qto_base <- qto_sing
    pre_base <- pre_sing
    sw1_base <- sw1_sing
    sw2_base <- sw2_sing
    sw3_base <- sw3_sing
  }else{
    epn_base <- cbind(epn_base, epn_sing)
    sno_base <- cbind(sno_base, sno_sing)
    qto_base <- cbind(qto_base, qto_sing)
    pre_base <- cbind(pre_base, pre_sing)
    sw1_base <- cbind(sw1_base, sw1_sing)
    sw2_base <- cbind(sw2_base, sw2_sing)
    sw3_base <- cbind(sw3_base, sw3_sing)
  }
  
}

#get simulation results for Cochem watershed
for (i in 1:length(cube_index_col_coch)) {
  
  print(paste(i, "of", length(cube_index_col_coch)))
  
  epn_sing <- pef_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  sno_sing <- snow_cube[cube_index_col_coch[i], cube_index_row_coch[i], ]
  qto_sing <- qto_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  pre_sing <- pre_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  sw1_sing <- sw1_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  sw2_sing <- sw2_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  sw3_sing <- sw3_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
  
  if(i == 1){
    epn_coch <- epn_sing
    sno_coch <- sno_sing
    qto_coch <- qto_sing
    pre_coch <- pre_sing
    sw1_coch <- sw1_sing
    sw2_coch <- sw2_sing
    sw3_coch <- sw3_sing
  }else{
    epn_coch <- cbind(epn_coch, epn_sing)
    sno_coch <- cbind(sno_coch, sno_sing)
    qto_coch <- cbind(qto_coch, qto_sing)
    pre_coch <- cbind(pre_coch, pre_sing)
    sw1_coch <- cbind(sw1_coch, sw1_sing)
    sw2_coch <- cbind(sw2_coch, sw2_sing)
    sw3_coch <- cbind(sw3_coch, sw3_sing)
  }
  
}

#Values on basin scale
base_ep_mea <- c(NA, apply(epn_base, 1, mea_na)) #fluxes/states only start from 02.01.1954
base_qt_mea <- c(NA, apply(qto_base, 1, mea_na))
base_sd_mea <- c(NA, apply(sno_base, 1, mea_na))
base_w1_mea <- c(NA, apply(sw1_base, 1, mea_na))
base_w2_mea <- c(NA, apply(sw2_base, 1, mea_na))
base_w3_mea <- c(NA, apply(sw3_base, 1, mea_na))
base_pr_mea <- apply(pre_base, 1, mea_na)
base_sd_mea_dif <- c(NA, diff(base_sd_mea))
base_sd_mea_dif[which(base_sd_mea_dif > 0)] <- NA
base_me_mea <- base_sd_mea_dif * -1 #melt positive values
base_wc_mea <- apply(cbind(base_w1_mea, base_w2_mea, base_w3_mea), 1, mea_na)

coch_ep_mea <- c(NA, apply(epn_coch, 1, mea_na)) #fluxes/states only start from 02.01.1954
coch_qt_mea <- c(NA, apply(qto_coch, 1, mea_na))
coch_sd_mea <- c(NA, apply(sno_coch, 1, mea_na))
coch_w1_mea <- c(NA, apply(sw1_coch, 1, mea_na))
coch_w2_mea <- c(NA, apply(sw2_coch, 1, mea_na))
coch_w3_mea <- c(NA, apply(sw3_coch, 1, mea_na))
coch_pr_mea <- apply(pre_coch, 1, mea_na)
coch_sd_mea_dif <- c(NA, diff(coch_sd_mea))
coch_sd_mea_dif[which(coch_sd_mea_dif > 0)] <- NA
coch_me_mea <- coch_sd_mea_dif * -1 #melt positive values
coch_wc_mea <- apply(cbind(coch_w1_mea, coch_w2_mea, coch_w3_mea), 1, mea_na)

#Mean annual cycle volumetric soil water content first 30 cm
base_wc_mea_day <- ord_day(base_wc_mea,
                           date = date_simu,
                           start_y = 1954,
                           end_y = 2013)
base_wc_ann <- apply(base_wc_mea_day, 2, mea_na)

coch_wc_mea_day <- ord_day(coch_wc_mea,
                           date = date_simu,
                           start_y = 1954,
                           end_y = 2013)
coch_wc_ann <- apply(coch_wc_mea_day, 2, mea_na)

#Moving window snowmelt
base_melt_ma <- rollapply(data = base_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
coch_melt_ma <- rollapply(data = coch_me_mea, width = 14,
                          FUN = sum_na, align = "right", fill = NA)
#Moving window liquid precipitation
base_prec <- base_ep_mea - base_me_mea
coch_prec <- coch_ep_mea - coch_me_mea
base_prec_ma <- rollapply(data = base_prec, width = 5,
                             FUN = sum_na, align = "right", fill = NA)
coch_prec_ma <- rollapply(data = coch_prec, width = 5,
                          FUN = sum_na, align = "right", fill = NA)

#Moving window total precipitation
base_pret_ma <- rollapply(data = base_pr_mea, width = 5,
                          FUN = sum_na, align = "right", fill = NA)
coch_pret_ma <- rollapply(data = coch_pr_mea, width = 5,
                          FUN = sum_na, align = "right", fill = NA)

#Get runoff peaks
quan_thres <- 0.90
number_peaks <- 60
time_cond <- 21

thres_val_base <- quantile(simu_base, quan_thres, na.rm = T)
thres_val_coch <- quantile(simu_coch, quan_thres, na.rm = T)
thres_val_base_melt <- quantile(base_melt_ma, quan_thres, na.rm = T)
thres_val_coch_melt <- quantile(coch_melt_ma, quan_thres, na.rm = T)
thres_val_base_prli <- quantile(base_prec_ma, quan_thres, na.rm = T)
thres_val_coch_prli <- quantile(coch_prec_ma, quan_thres, na.rm = T)
thres_val_base_prto <- quantile(base_pret_ma, quan_thres, na.rm = T)
thres_val_coch_prto <- quantile(coch_pret_ma, quan_thres, na.rm = T)

pot_data_base <- data.frame(obs  = simu_base,
                            time = date_simu)
pot_data_coch <- data.frame(obs  = simu_coch,
                            time = date_simu)
pot_data_base_melt <- data.frame(obs  = base_melt_ma,
                            time = date_simu)
pot_data_coch_melt <- data.frame(obs  = coch_melt_ma,
                            time = date_simu)
pot_data_base_prli <- data.frame(obs  = base_prec_ma,
                                 time = date_simu)
pot_data_coch_prli <- data.frame(obs  = coch_prec_ma,
                                 time = date_simu)
pot_data_base_prto <- data.frame(obs  = base_pret_ma,
                                 time = date_simu)
pot_data_coch_prto <- data.frame(obs  = coch_pret_ma,
                                 time = date_simu)

pot_peaks_base_all <- clust(data = pot_data_base, u = thres_val_base, 
                        tim.cond = time_cond, clust.max = T, plot = F)
pot_peaks_coch_all <- clust(data = pot_data_coch, u = thres_val_coch, 
                            tim.cond = time_cond, clust.max = T, plot = F)
pot_data_base_melt$obs[which(is.na(pot_data_base_melt$obs))] <- 0
pot_data_coch_melt$obs[which(is.na(pot_data_coch_melt$obs))] <- 0
pot_peaks_base_all_melt <- clust(data = pot_data_base_melt, u = thres_val_base_melt, 
                            tim.cond = time_cond, clust.max = T, plot = F)
pot_peaks_coch_all_melt <- clust(data = pot_data_coch_melt, u = thres_val_coch_melt, 
                            tim.cond = time_cond, clust.max = T, plot = F)
pot_data_base_prli$obs[which(is.na(pot_data_base_prli$obs))] <- 0
pot_data_coch_prli$obs[which(is.na(pot_data_coch_prli$obs))] <- 0
pot_peaks_base_all_prli <- clust(data = pot_data_base_prli, u = thres_val_base_prli, 
                                 tim.cond = time_cond, clust.max = T, plot = F)
pot_peaks_coch_all_prli <- clust(data = pot_data_coch_prli, u = thres_val_coch_prli, 
                                 tim.cond = time_cond, clust.max = T, plot = F)
pot_data_base_prto$obs[which(is.na(pot_data_base_prto$obs))] <- 0
pot_data_coch_prto$obs[which(is.na(pot_data_coch_prto$obs))] <- 0
pot_peaks_base_all_prto <- clust(data = pot_data_base_prto, u = thres_val_base_prto, 
                                 tim.cond = time_cond, clust.max = T, plot = F)
pot_peaks_coch_all_prto <- clust(data = pot_data_coch_prto, u = thres_val_coch_prto, 
                                 tim.cond = time_cond, clust.max = T, plot = F)

pot_peaks_base_all_ord <- pot_peaks_base_all[order(pot_peaks_base_all[, 2], decreasing = T), ]
pot_peaks_coch_all_ord <- pot_peaks_coch_all[order(pot_peaks_coch_all[, 2], decreasing = T), ]
pot_peaks_base_all_melt_ord <- pot_peaks_base_all_melt[order(pot_peaks_base_all_melt[, 2], decreasing = T), ]
pot_peaks_coch_all_melt_ord <- pot_peaks_coch_all_melt[order(pot_peaks_coch_all_melt[, 2], decreasing = T), ]
pot_peaks_base_all_prli_ord <- pot_peaks_base_all_prli[order(pot_peaks_base_all_prli[, 2], decreasing = T), ]
pot_peaks_coch_all_prli_ord <- pot_peaks_coch_all_prli[order(pot_peaks_coch_all_prli[, 2], decreasing = T), ]
pot_peaks_base_all_prto_ord <- pot_peaks_base_all_prto[order(pot_peaks_base_all_prto[, 2], decreasing = T), ]
pot_peaks_coch_all_prto_ord <- pot_peaks_coch_all_prto[order(pot_peaks_coch_all_prto[, 2], decreasing = T), ]

pot_peaks_base <- pot_peaks_base_all_ord[1:number_peaks, ]
pot_peaks_coch <- pot_peaks_coch_all_ord[1:number_peaks, ]
pot_peaks_base_melt <- pot_peaks_base_all_melt_ord[1:number_peaks, ]
pot_peaks_coch_melt <- pot_peaks_coch_all_melt_ord[1:number_peaks, ]
pot_peaks_base_prli <- pot_peaks_base_all_prli_ord[1:number_peaks, ]
pot_peaks_coch_prli <- pot_peaks_coch_all_prli_ord[1:number_peaks, ]
pot_peaks_base_prto <- pot_peaks_base_all_prto_ord[1:number_peaks, ]
pot_peaks_coch_prto <- pot_peaks_coch_all_prto_ord[1:number_peaks, ]

melt_peak_base <- base_melt_ma[pot_peaks_base[, 3]]
prec_peak_base <- base_prec_ma[pot_peaks_base[, 3]]
melt_peak_coch <- coch_melt_ma[pot_peaks_coch[, 3]]
prec_peak_coch <- coch_prec_ma[pot_peaks_coch[, 3]]

#fraction snowmelt contribution
peak_frac_mel_base <- melt_peak_base / (melt_peak_base + prec_peak_base)
peak_frac_mel_coch <- melt_peak_coch / (melt_peak_coch + prec_peak_coch)

#day of the year of events
peak_doy_base <- as.numeric(format(date_simu[pot_peaks_base[, 3]], '%j'))
peak_doy_coch <- as.numeric(format(date_simu[pot_peaks_coch[, 3]], '%j'))
peak_doy_base_melt <- as.numeric(format(date_simu[pot_peaks_base_melt[, 3]], '%j'))
peak_doy_coch_melt <- as.numeric(format(date_simu[pot_peaks_coch_melt[, 3]], '%j'))
peak_doy_base_prli <- as.numeric(format(date_simu[pot_peaks_base_prli[, 3]], '%j'))
peak_doy_coch_prli <- as.numeric(format(date_simu[pot_peaks_coch_prli[, 3]], '%j'))
peak_doy_base_prto <- as.numeric(format(date_simu[pot_peaks_base_prto[, 3]], '%j'))
peak_doy_coch_prto <- as.numeric(format(date_simu[pot_peaks_coch_prto[, 3]], '%j'))

mel_thres_1 <- 0.33
mel_thres_2 <- 0.66
mel_1_ind_base <- which(peak_frac_mel_base > mel_thres_1)
mel_2_ind_base <- which(peak_frac_mel_base > mel_thres_2)
mel_1_ind_coch <- which(peak_frac_mel_coch > mel_thres_1)
mel_2_ind_coch <- which(peak_frac_mel_coch > mel_thres_2)



pdf(paste0(bas_dir, "res_figs/runoff_seas",".pdf"), width = 12.0, height = 12)

par(family = "serif")

layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
              7, 2, byrow = T), widths=c(), heights=c(0.2, 1, 1, 1, 1, 1, 1))

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

cex_points <- 1.9
cex_main <- 1.3
cex_x_label <- 1.5
alpha_sel <- 0.5

par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("Basel",  side = 3, line = -2.2, cex = 1.5, adj = 0.5, outer = F)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("Cochem",  side = 3, line = -2.2, cex = 1.5, adj = 0.5, outer = F)

par(mar = c(2, 3.5, 2.5, 0.5))

#Basel runoff peaks
plot(peak_doy_base, pot_peaks_base[, 2], pch = 19, ylab = "", xlab = "", axes = F,
     xlim = c(0, 365), type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_base, pot_peaks_base[, 2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("a) Runoff peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Cochem runoff peaks
plot(peak_doy_coch, pot_peaks_coch[, 2], pch = 19, ylab = "", xlab = "", axes = F,
     xlim = c(0, 365), type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_coch, pot_peaks_coch[, 2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("b) Runoff peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Basel snowmelt fraction
plot(peak_doy_base, peak_frac_mel_base, pch = 19, ylab = "", xlab = "", axes = F,
     xlim = c(0, 365), ylim = c(0, 1), type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_base, peak_frac_mel_base, col = alpha("black", alpha = alpha_sel), 
       pch = 19, cex = cex_points)
# points(peak_doy_base[mel_1_ind_base], peak_frac_mel_base[mel_1_ind_base], col = alpha("black", alpha = alpha_sel), 
#        pch = 19, cex = cex_points)
# points(peak_doy_base[mel_2_ind_base], peak_frac_mel_base[mel_2_ind_base], col = alpha("black", alpha = alpha_sel), 
#        pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("c) Snowmelt contribution", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("Fraction melt [-]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Cochem snowmelt fraction
plot(peak_doy_coch, peak_frac_mel_coch, pch = 19, ylab = "", xlab = "", axes = F,
     xlim = c(0, 365), ylim = c(0, 1), type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_coch, peak_frac_mel_coch, col = alpha("black", alpha = alpha_sel), 
       pch = 19, cex = cex_points)
# points(peak_doy_coch[mel_1_ind_coch], peak_frac_mel_coch[mel_1_ind_coch], col = alpha("black", alpha = alpha_sel), 
#        pch = 19, cex = cex_points)
# points(peak_doy_coch[mel_2_ind_coch], peak_frac_mel_coch[mel_2_ind_coch], col = alpha("black", alpha = alpha_sel), 
#        pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("d) Snowmelt contribution", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("Fraction melt [-]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")


#Basel melt peaks 
plot(peak_doy_base_melt, pot_peaks_base_melt[ ,2], type = "n", xlim = c(0, 365), ylab = "", xlab = "", axes = F, 
     ylim = range(c(pot_peaks_base_melt[ ,2], pot_peaks_coch_melt[ ,2])))
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_base_melt, pot_peaks_base_melt[ , 2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("e) Snowmelt peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("14-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")


#Cochem melt peaks
plot(peak_doy_coch_melt, pot_peaks_coch_melt[ ,2], type = "n", xlim = c(0, 365), ylab = "", xlab = "", axes = F,
     ylim = range(c(pot_peaks_base_melt[ ,2], pot_peaks_coch_melt[ ,2])))
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_coch_melt, pot_peaks_coch_melt[ ,2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("f) Snowmelt peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("14-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")


#Basel liquid precipitation
plot(peak_doy_base_prli, pot_peaks_base_prli[ ,2], xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_base_prli, pot_peaks_base_prli[ ,2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("g) Liquid rainfall peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("5-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Cochem liquid precipiation
plot(peak_doy_coch_prli, pot_peaks_coch_prli[ ,2], xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_coch_prli, pot_peaks_coch_prli[ ,2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("h) Liquid rainfall peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("5-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Basel total precipitation
plot(peak_doy_base_prto, pot_peaks_base_prto[ ,2], xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_base_prto, pot_peaks_base_prto[ ,2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("i) Total rainfall peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("5-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Cochem total precipitation
plot(peak_doy_coch_prto, pot_peaks_coch_prto[ ,2], xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n")
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
points(peak_doy_coch_prto, pot_peaks_coch_prto[ ,2], col = alpha("black", alpha = alpha_sel), pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("j) Total rainfall peaks", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("5-day sum [mm]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Basel Soil water content
plot(base_wc_ann, xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n",
     ylim = c(min_na(c(base_wc_ann, coch_wc_ann)), 1))
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
lines(base_wc_ann, col = "black", pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("k) Soil water content (0-30 cm)", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("Vol. SWC [-]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

#Cochem Soil water content
plot(coch_wc_ann, xlim = c(0, 365), ylab = "", xlab = "", axes = F, type = "n", 
     ylim = c(min_na(c(base_wc_ann, coch_wc_ann)), 1))
abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
lines(coch_wc_ann, col = "black", pch = 19, cex = cex_points)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = cex_x_label)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
mtext("l) Soil water content (0-30 cm)", side = 3, line = 0.3, cex = cex_main, adj = 0.0)
mtext("Vol. SWC [-]", side = 2, line = 1.6, cex = 1.1)
box(bty = "o")

dev.off()





#Runoff timing
vio_dat <- data.frame(doy = c(peak_doy_coch, peak_doy_base),
                      basin = c(rep("Cochem", length(peak_doy_coch)), 
                                rep("Basel",  length(peak_doy_base))))
vio_dat$basin <- as.factor(vio_dat$basin)

vio_run <- ggplot(vio_dat, aes(x = basin, y = doy, fill = basin)) +
  geom_violin(trim = F, show.legend = F) +
  ylim(0, 365) +
  coord_flip() +
  scale_fill_manual(values=c("grey40", "grey75")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        text = element_text(size = 12)) +
  labs(title="a) Runoff timing", y="Day of the year", x = "") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 7, 
               show.legend = F, fill = "black")

#Runoff histograms
his_dat <- data.frame(mag = c(pot_peaks_coch[, 2], pot_peaks_base[, 2]),
                      basin = c(rep("Cochem", length(peak_doy_coch)), 
                                rep("Basel",  length(peak_doy_base))))
his_dat$basin <- as.factor(his_dat$basin)

his_run <- ggplot(his_dat, aes(x = mag, fill = basin)) +
  geom_histogram(alpha = 1.0, binwidth = 250) +
  scale_fill_manual(values=c("grey40", "grey75")) +
  labs(title="b) Runff magnitudes", y="", x = "Magnitude [mm/14d]") +
  theme_bw() +
  theme(text = element_text(size = 12))

#Snowmelt timing
vio_dat <- data.frame(doy = c(peak_doy_coch_melt, peak_doy_base_melt),
                      basin = c(rep("Cochem", length(peak_doy_coch_melt)), 
                                rep("Basel",  length(peak_doy_base_melt))))
vio_dat$basin <- as.factor(vio_dat$basin)

vio_sno <- ggplot(vio_dat, aes(x = basin, y = doy, fill = basin)) +
  geom_violin(trim = F, show.legend = F) +
  ylim(0, 365) +
  coord_flip() +
  scale_fill_manual(values=c("grey40", "grey75")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        text = element_text(size = 12)) +
  labs(title="c) 14d snowmelt timing", y="Day of the year", x = "") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 7, 
               show.legend = F, fill = "black")


#Snowmelt histograms
his_dat <- data.frame(mag = c(pot_peaks_coch_melt[, 2], pot_peaks_base_melt[, 2]),
                      basin = c(rep("Cochem", length(peak_doy_coch_melt)), 
                                rep("Basel",  length(peak_doy_base_melt))))
his_dat$basin <- as.factor(his_dat$basin)

his_sno <- ggplot(his_dat, aes(x = mag, fill = basin)) +
  geom_histogram(alpha = 1.0, binwidth = 5) +
  scale_fill_manual(values=c("grey40", "grey75")) +
  labs(title="d) 14d snowmelt magnitudes", y="", x = "Magnitude [mm/14d]") +
  theme_bw() +
  theme(text = element_text(size = 12))

#Liquid rainfall timing
vio_dat <- data.frame(doy = c(peak_doy_coch_prli, peak_doy_base_prli),
                      basin = c(rep("Cochem", length(peak_doy_coch_prli)), 
                                rep("Basel",  length(peak_doy_base_prli))))
vio_dat$basin <- as.factor(vio_dat$basin)

vio_rli <- ggplot(vio_dat, aes(x = basin, y = doy, fill = basin)) +
  geom_violin(trim = F, show.legend = F) +
  ylim(0, 365) +
  coord_flip() +
  scale_fill_manual(values=c("grey40", "grey75")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        text = element_text(size = 12)) +
  labs(title="e) 5d liquid rainfall timing", y="Day of the year", x = "") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 7, 
               show.legend = F, fill = "black")


#Liquid rainfall histograms
his_dat <- data.frame(mag = c(pot_peaks_coch_prli[, 2], pot_peaks_base_prli[, 2]),
                      basin = c(rep("Cochem", length(peak_doy_coch_prli)), 
                                rep("Basel",  length(peak_doy_base_prli))))
his_dat$basin <- as.factor(his_dat$basin)

his_rli <- ggplot(his_dat, aes(x = mag, fill = basin)) +
  geom_histogram(alpha = 1.0, binwidth = 5) +
  scale_fill_manual(values=c("grey40", "grey75")) +
  labs(title="f) 5d liquid rainfall magnitudes", y="", x = "Magnitude [mm/5d]") +
  theme_bw() +
  theme(text = element_text(size = 12))

#Total rainfall timing
vio_dat <- data.frame(doy = c(peak_doy_coch_prto, peak_doy_base_prto),
                      basin = c(rep("Cochem", length(peak_doy_coch_prto)), 
                                rep("Basel",  length(peak_doy_base_prto))))
vio_dat$basin <- as.factor(vio_dat$basin)

vio_rto <- ggplot(vio_dat, aes(x = basin, y = doy, fill = basin)) +
  geom_violin(trim = F, show.legend = F) +
  ylim(0, 365) +
  coord_flip() +
  scale_fill_manual(values=c("grey40", "grey75")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        text = element_text(size = 12)) +
  labs(title="g) 5d total rainfall timing", y="Day of the year", x = "") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 7, 
               show.legend = F, fill = "black")


#Liquid rainfall histograms
his_dat <- data.frame(mag = c(pot_peaks_coch_prto[, 2], pot_peaks_base_prto[, 2]),
                      basin = c(rep("Cochem", length(peak_doy_coch_prto)), 
                                rep("Basel",  length(peak_doy_base_prto))))
his_dat$basin <- as.factor(his_dat$basin)

his_rto <- ggplot(his_dat, aes(x = mag, fill = basin)) +
  geom_histogram(alpha = 1.0, binwidth = 5) +
  scale_fill_manual(values=c("grey40", "grey75")) +
  labs(title="h) 5d total rainfall magnitudes", y="", x = "Magnitude [mm/5d]") +
  theme_bw() +
  theme(text = element_text(size = 12))

#Snowmelt contribution timing
his_dat <- data.frame(mag = c(peak_frac_mel_coch, peak_frac_mel_base),
                      doy = c(peak_doy_coch, peak_doy_base),
                      basin = c(rep("Cochem", length(peak_frac_mel_coch)), 
                                rep("Basel",  length(peak_frac_mel_base))))
his_dat$basin <- as.factor(his_dat$basin)

his_cot <- ggplot(his_dat, aes(x = doy, y = mag, color = basin)) +
  geom_point(aes(size=mag), show.legend = T) +
  scale_color_manual(values=c("grey40", "grey75")) +
  labs(title="i) Snowmelt contribution timing", x = "Day of the year", y = "Fraction [-]") +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 12))

#Snowmelt contribution histograms
his_com <- ggplot(his_dat, aes(x = mag, fill = basin)) +
  geom_histogram(alpha = 1.0, binwidth = 0.1, na.rm = T) +
  scale_fill_manual(values=c("grey40", "grey75")) +
  labs(title="j) Snowmelt contribution magnitudes", y="", x = "Fraction [-]") +
  theme_bw() +
  theme(text = element_text(size = 12))

#Soil water content Cochem
swc_dat <- data.frame(mag = c(coch_wc_mea),
                      doy = as.numeric(format(date_simu, "%j")))
swc_dat$doy <- as.factor(swc_dat$doy)

swc_dat_peak <- data.frame(mag = coch_wc_mea[pot_peaks_coch[, 3]],
                           doy = peak_doy_coch)
swc_dat_peak$doy <- as.factor(swc_dat_peak$doy)

swc_coch <- ggplot(swc_dat, aes(x = doy, y = mag)) +
  geom_point(alpha = 0.2, size = 1, shape = 16) +
  geom_point(data = swc_dat_peak, aes(x = doy, y = mag), colour = "red", alpha = 1.0) +
  labs(title="k) Cochem soil water content (0-30 cm)", x = "Day of the year", y = "Vol. SWC [-]") +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 12)) +
  scale_x_discrete(breaks=seq(50, 350, 50), labels=seq(50, 350, 50)) +
  ylim(min_na(c(base_wc_mea, coch_wc_mea)), 1)

#Soil water content Basel
swc_dat <- data.frame(mag = c(base_wc_mea),
                      doy = as.numeric(format(date_simu, "%j")))
swc_dat$doy <- as.factor(swc_dat$doy)

swc_dat_peak <- data.frame(mag = base_wc_mea[pot_peaks_base[, 3]],
                           doy = peak_doy_base)
swc_dat_peak$doy <- as.factor(swc_dat_peak$doy)

swc_base <- ggplot(swc_dat, aes(x = doy, y = mag)) +
  geom_point(alpha = 0.2, size = 1, shape = 16) +
  geom_point(data = swc_dat_peak, aes(x = doy, y = mag), colour = "red", alpha = 1.0) +
  labs(title="l) Basel soil water content (0-30 cm)", x = "Day of the year", y = "Vol. SWC [-]") +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 12)) +
  scale_x_discrete(breaks=seq(50, 350, 50), labels=seq(50, 350, 50)) +
  ylim(min_na(c(base_wc_mea, coch_wc_mea)), 1)


pdf(paste0(bas_dir, "res_figs/runoff_seaso",".pdf"), width = 12.0, height = 15)

ggarrange(vio_run, his_run, 
          vio_sno, his_sno,
          vio_rli, his_rli,
          vio_rto, his_rto,
          his_cot, his_com,
          swc_coch, swc_base,
          ncol = 2, nrow = 6,
          heights = c(1, 1))

dev.off()


#seas_flod_old----

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
qto_mea <- apply(qto_cube, c(1, 2), sum_na) / length(sta_yea:end_yea) #[mm]

snow_max_c <- c(snow_max)
pef_mea_c <- c(pef_mea)
qto_mea_c <- c(qto_mea)

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

cols_spat_qto <- foreach(i = 1:length(qto_mea_c), .combine = 'cbind') %dopar% {
  
  val2col(val_in = qto_mea_c[i],
          dat_ref = qto_mea_c,
          do_bicol = F,
          virid_dir = -1,
          do_log = F,
          col_na = "white")
  
}

#Plot times series with POT
pdf(paste0(bas_dir,"res_figs/map_flood.pdf"), width = 16, height = 6)

layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4,  rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())
# layout.show(n=6)

par(family = "serif")
cex_pch <- 1.20
mar_1 <- c(1.5, 0.5, 1.5, 0.5)

#Plot: Snowpack
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_sno, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("a) Snow depth", side = 3, line = -1.0, cex = 1.5)

rect(xleft = 7.1, xright = 10.3, ybottom = 46.2, ytop = 47.5)
rect(xleft = 5.3, xright = 12.0, ybottom = 48.1, ytop = 50.7)
text(10.5, 50.8, "pluvial", cex = 1.8)
text(10.5, 47.6, "nival", cex = 1.8)

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
mtext("b) Effective Precipitation", side = 3, line = -1.0, cex = 1.5)

rect(xleft = 7.1, xright = 10.3, ybottom = 46.2, ytop = 47.5)
rect(xleft = 5.3, xright = 12.0, ybottom = 48.1, ytop = 50.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(pef_mea_c, na.rm = T)[1], range(pef_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(pef_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

#Plot: Total discharge generated per cell
par(mar = mar_1)
plot(c(lon), c(lat), pch = 15, col = cols_spat_qto, cex = 1.0, axes = F, ylab = "", xlab = "")
mtext("c) Discharge generated", side = 3, line = -1.0, cex = 1.5)

rect(xleft = 7.1, xright = 10.3, ybottom = 46.2, ytop = 47.5)
rect(xleft = 5.3, xright = 12.0, ybottom = 48.1, ytop = 50.7)

par(mar = c(2.0, 0.2, 5.0, 2.9))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(range(qto_mea_c, na.rm = T)[1], range(qto_mea_c, na.rm = T)[2], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(qto_mea_c), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext("[mm]", side = 3, line = 0.7, cex = 1.2)
box()

dev.off()




#Select cells
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
  sno_sing <- snow_cube[cube_index_col_2[i], cube_index_row_2[i], ]
  qto_sing <- qto_cube [cube_index_col_2[i], cube_index_row_2[i], ]
  
  if(i == 1){
    epns_plu <- epn_sing
    snos_plu <- sno_sing
    qtos_plu <- qto_sing
  }else{
    epns_plu <- cbind(epns_plu, epn_sing)
    snos_plu <- cbind(snos_plu, sno_sing)
    qtos_plu <- cbind(qtos_plu, qto_sing)
  }
}

plu_ep_sum <- apply(epns_plu, 1, sum_na)
plu_qt_sum <- apply(qtos_plu, 1, sum_na)
plu_sd_sum <- apply(snos_plu, 1, sum_na)
plu_sd_sum_dif <- c(NA, diff(plu_sd_sum))
plu_sd_sum_dif[which(plu_sd_sum_dif > 0)] <- NA
plu_sn_sum <- plu_sd_sum_dif * -1 #melt positive values


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

window_plu_sn <- 5
plu_sn_sum_ma <- rollapply(data = plu_sn_sum, width = window_plu_sn,
                           FUN = sum_na, align = "center", fill = NA)

plot(plu_ep_sum_ma, type = "l")
plot(plu_qt_sum_ma, type = "l")
plot(plu_sn_sum_ma, type = "l")

#Peak over threshold
pot_thre_niv_sn <- quantile(niv_sn_sum_ma, 0.95, na.rm = T)
pot_thre_niv_ep <- quantile(niv_ep_sum_ma, 0.95, na.rm = T)
pot_thre_niv_qt <- quantile(niv_qt_sum_ma, 0.95, na.rm = T)
pot_thre_plu_ep <- quantile(plu_ep_sum_ma, 0.95, na.rm = T)
pot_thre_plu_qt <- quantile(plu_qt_sum_ma, 0.95, na.rm = T)
pot_thre_plu_sn <- quantile(plu_sn_sum_ma, 0.95, na.rm = T)

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
pot_data_plu_sn <- data.frame(obs = plu_sn_sum_ma,
                              time = date[-1])

pot_data_niv_sn$obs[is.na(pot_data_niv_sn$obs)] <- 0
pot_data_niv_ep$obs[is.na(pot_data_niv_ep$obs)] <- 0
pot_data_niv_qt$obs[is.na(pot_data_niv_qt$obs)] <- 0
pot_data_plu_ep$obs[is.na(pot_data_plu_ep$obs)] <- 0
pot_data_plu_qt$obs[is.na(pot_data_plu_qt$obs)] <- 0
pot_data_plu_sn$obs[is.na(pot_data_plu_sn$obs)] <- 0

pot_peaks_niv_sn <- clust(data = pot_data_niv_sn, u = pot_thre_niv_sn, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_niv_ep <- clust(data = pot_data_niv_ep, u = pot_thre_niv_ep, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_niv_qt <- clust(data = pot_data_niv_qt, u = pot_thre_niv_qt, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_plu_ep <- clust(data = pot_data_plu_ep, u = pot_thre_plu_ep, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_plu_qt <- clust(data = pot_data_plu_qt, u = pot_thre_plu_qt, tim.cond = 14, clust.max = T, plot = F)
pot_peaks_plu_sn <- clust(data = pot_data_plu_sn, u = pot_thre_plu_sn, tim.cond = 14, clust.max = T, plot = F)

#Plot times series with POT
pdf(paste0(bas_dir,"res_figs/pot_flood.pdf"), width = 16, height = 6)

mar_sea <- c(2, 5.0, 3, 0.5)
par(mar = mar_sea)
par(mfrow = c(2, 3))
par(family = "serif")

plot(pot_data_niv_sn$time, niv_sn_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_niv_sn$time[pot_peaks_niv_sn[, 3]], pot_peaks_niv_sn[, 2], col = "darkred", 
       cex = 1.2, pch = 19)
mtext("Moving 14d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("a) Nival area - Snowmelt", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

plot(pot_data_niv_ep$time, niv_ep_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_niv_ep$time[pot_peaks_niv_ep[, 3]], pot_peaks_niv_ep[, 2], col = "orange3",
       cex = 1.2, pch = 19)
mtext("Moving 14d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("b) Nival area - Effective precip.", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

plot(pot_data_niv_qt$time, niv_qt_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_niv_ep$time[pot_peaks_niv_qt[, 3]], pot_peaks_niv_qt[, 2], col = "gold3",
       cex = 1.2, pch = 19)
mtext("Moving 14d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("c) Nival area - Discharge generated", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

plot(pot_data_plu_sn$time, plu_sn_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_plu_sn$time[pot_peaks_plu_sn[, 3]], pot_peaks_plu_sn[, 2], col = "black",
       cex = 1.2, pch = 19)
mtext("Moving 5d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("d) Pluvial area - Snowmelt", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

plot(pot_data_plu_ep$time, plu_ep_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_plu_sn$time[pot_peaks_plu_ep[, 3]], pot_peaks_plu_ep[, 2], col = "darkblue",
       cex = 1.2, pch = 19)
mtext("Moving 5d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("e) Pluvial area - Effective precip.", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

plot(pot_data_plu_qt$time, plu_qt_sum_ma, type = "l", ylab = "", xlab = "", cex.axis = 1.2)
points(pot_data_plu_qt$time[pot_peaks_plu_qt[, 3]], pot_peaks_plu_qt[, 2], col = "steelblue3",
       cex = 1.2, pch = 19)
mtext("Moving 5d sum [mm]", side = 2, line = 2.4, cex = 1.2, adj = 0.5)
mtext("f) Pluvial area - Discharge generated", side = 3, line = 0.5, cex = 1.2, adj = 0.0)

dev.off()


peaks_doy_niv_sn <- as.numeric(format(pot_data_niv_sn$time[pot_peaks_niv_sn[, 3]], '%j'))
peaks_doy_niv_ep <- as.numeric(format(pot_data_niv_ep$time[pot_peaks_niv_ep[, 3]], '%j'))
peaks_doy_niv_qt <- as.numeric(format(pot_data_niv_qt$time[pot_peaks_niv_qt[, 3]], '%j'))
peaks_doy_plu_ep <- as.numeric(format(pot_data_plu_ep$time[pot_peaks_plu_ep[, 3]], '%j'))
peaks_doy_plu_qt <- as.numeric(format(pot_data_plu_qt$time[pot_peaks_plu_qt[, 3]], '%j'))
peaks_doy_plu_sn <- as.numeric(format(pot_data_plu_sn$time[pot_peaks_plu_sn[, 3]], '%j'))


pdf(paste0(bas_dir,"res_figs/seas_flood.pdf"), width = 16, height = 6)

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

ylims <- range(c(pot_peaks_niv_sn[ ,2], pot_peaks_niv_ep[ ,2], pot_peaks_niv_qt[ ,2], 
                 pot_peaks_plu_ep[, 2], pot_peaks_plu_qt[, 2], pot_peaks_plu_sn[, 2]))

mar_sea <- c(2, 3.5, 3, 0.5)
par(mar = mar_sea)
par(mfrow = c(2, 3))
par(family = "serif")
alpha_sel <- 0.4

#Nival part - Snowmelt
plot(peaks_doy_niv_sn, pot_peaks_niv_sn[, 2], col = alpha("darkred", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("a) Nival area - Snowmelt", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

#Nival part - Effective precipitation
plot(peaks_doy_niv_ep, pot_peaks_niv_ep[, 2], col = alpha("orange3", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("a) Nival area - Effective precip.", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

#Nival part - Discharge generated
plot(peaks_doy_niv_qt, pot_peaks_niv_qt[, 2], col = alpha("gold3", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("c) Nival area - Discharge generated", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

#Pluvial part - Snowmelt
plot(peaks_doy_plu_sn, pot_peaks_plu_sn[, 2], col = alpha("black", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("d) Pluvial area - Snowmelt", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

#Pluvial part - Effective precipitation
plot(peaks_doy_plu_ep, pot_peaks_plu_ep[, 2], col = alpha("darkblue", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("e) Pluvial area - Effective precip.", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

#Pluvial part - Discharge generated
plot(peaks_doy_plu_qt, pot_peaks_plu_qt[, 2], col = alpha("steelblue3", alpha = alpha_sel), pch = 19, 
     cex = 2, xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = ylims)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.2)
mtext("Event magnitude [mm]", side = 2, line = 1.7, cex = 1.2, adj = 0.5)
mtext("f) Pluvial area - Discharge generated", side = 3, line = 0.5, cex = 1.2, adj = 0.0)
box()

dev.off()



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

#seas_plot----

grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

lobi_file <- paste0(grdc_dir, "6435060_Q_Day.Cmd.txt")
koel_file <- paste0(grdc_dir, "6335060_Q_Day.Cmd.txt")
coch_file <- paste0(grdc_dir, "6336050_Q_Day.Cmd.txt")
kaub_file <- paste0(grdc_dir, "6335100_Q_Day.Cmd.txt")
wuer_file <- paste0(grdc_dir, "6335500_Q_Day.Cmd.txt")
worm_file <- paste0(grdc_dir, "6335180_Q_Day.Cmd.txt")
rock_file <- paste0(grdc_dir, "6335600_Q_Day.Cmd.txt")
spey_file <- paste0(grdc_dir, "6335170_Q_Day.Cmd.txt")
base_file <- paste0(grdc_dir, "6935051_Q_Day.Cmd.txt")
unte_file <- paste0(grdc_dir, "6935300_Q_Day.Cmd.txt")
reki_file <- paste0(grdc_dir, "6935054_Q_Day.Cmd.txt")

disc_reki_full <- read_grdc(reki_file)
disc_unte_full <- read_grdc(unte_file)
disc_base_full <- read_grdc(base_file)
disc_spey_full <- read_grdc(spey_file)
disc_rock_full <- read_grdc(rock_file)
disc_worm_full <- read_grdc(worm_file)
disc_wuer_full <- read_grdc(wuer_file)
disc_kaub_full <- read_grdc(kaub_file)
disc_coch_full <- read_grdc(coch_file)
disc_koel_full <- read_grdc(koel_file)
disc_lobi_full <- read_grdc(lobi_file)

date_simu <- seq(as.Date("1954-01-01", format = "%Y-%m-%d"), 
                 as.Date("2013-12-31", format = "%Y-%m-%d"), by = "day")

disc_reki <- disc_reki_full[which(disc_reki_full$date %in% date_simu), ]
disc_unte <- disc_unte_full[which(disc_unte_full$date %in% date_simu), ]
disc_base <- disc_base_full[which(disc_base_full$date %in% date_simu), ]
disc_spey <- disc_spey_full[which(disc_spey_full$date %in% date_simu), ]
disc_rock <- disc_rock_full[which(disc_rock_full$date %in% date_simu), ]
disc_worm <- disc_worm_full[which(disc_worm_full$date %in% date_simu), ]
disc_wuer <- disc_wuer_full[which(disc_wuer_full$date %in% date_simu), ]
disc_kaub <- disc_kaub_full[which(disc_kaub_full$date %in% date_simu), ]
disc_coch <- disc_coch_full[which(disc_coch_full$date %in% date_simu), ]
disc_koel <- disc_koel_full[which(disc_koel_full$date %in% date_simu), ]
disc_lobi <- disc_lobi_full[which(disc_lobi_full$date %in% date_simu), ]

pdf(paste0(bas_dir, "res_figs/seas_plot",".pdf"), width = 5.0, height = 9.0)
# tiff(paste0(bas_dir, "res_figs/runoff_qu_", stat_sel,".tiff"), width = 3*2.0, height = 3*1.5,
#      units = "in", res = 800)
par(family = "serif")

par(mfrow = c(5, 1))
par(mar = c(2.5, 4.5, 3.0, 0.5))
width_sel <- 7

seas.var.plot(disc_reki, var = "value", width = width_sel, main = "a) Rekingen",
              ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))
seas.var.plot(disc_base, var = "value", width = width_sel, main = "b) Basel",
              ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))
seas.var.plot(disc_rock, var = "value", width = width_sel, main = "c) Rockenau",
              ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))
seas.var.plot(disc_coch, var = "value", width = width_sel, main = "d) Cochem",
              ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))
seas.var.plot(disc_koel, var = "value", width = width_sel, main = "e) Cologne",
              ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))

dev.off()


#coin_calc----

#Function binarize discharge data
bina_pot <- function(val_in, time_in, quan_thres = 0.90, do_numb = T, numb_events = 60, inde_con = 21){

  thres_val <- quantile(val_in, quan_thres, na.rm = T)
  
  pot_data <- data.frame(obs  = val_in,
                         time = time_in)
  
  pot_peaks <- clust(data = pot_data, u = thres_val, tim.cond = 14, clust.max = T, plot = F)
  
  if(do_numb){
  
    pot_peaks <- pot_peaks[order(pot_peaks[, 2], decreasing = T), ]
    pot_peaks <- pot_peaks[1:numb_events, ]
    
  }
  
  bina_vals <- rep(0, length(val_in))
  bina_vals[as.numeric(pot_peaks[, 3])] <- 1
  
  return(bina_vals)
  
}

bina_base_obs <- bina_pot(val_in = disc_base$value, time_in = disc_base$date)  
bina_spey_obs <- bina_pot(val_in = disc_spey$value, time_in = disc_spey$date)  
bina_rock_obs <- bina_pot(val_in = disc_rock$value, time_in = disc_rock$date)  
bina_worm_obs <- bina_pot(val_in = disc_worm$value, time_in = disc_worm$date)  
bina_wuer_obs <- bina_pot(val_in = disc_wuer$value, time_in = disc_wuer$date)  
bina_kaub_obs <- bina_pot(val_in = disc_kaub$value, time_in = disc_kaub$date)  
bina_coch_obs <- bina_pot(val_in = disc_coch$value, time_in = disc_coch$date)  
bina_koel_obs <- bina_pot(val_in = disc_koel$value, time_in = disc_koel$date)  
bina_lobi_obs <- bina_pot(val_in = disc_lobi$value, time_in = disc_lobi$date)  

bina_base_sim <- bina_pot(simu_base, date_simu)  
bina_spey_sim <- bina_pot(simu_spey, date_simu)  
bina_rock_sim <- bina_pot(simu_rock, date_simu)  
bina_worm_sim <- bina_pot(simu_worm, date_simu)  
bina_wuer_sim <- bina_pot(simu_wuer, date_simu)  
bina_kaub_sim <- bina_pot(simu_kaub, date_simu)  
bina_coch_sim <- bina_pot(simu_coch, date_simu)  
bina_koel_sim <- bina_pot(simu_koel, date_simu)  
bina_lobi_sim <- bina_pot(simu_lobi, date_simu)  

bina_obs <- cbind(bina_base_obs, bina_spey_obs, bina_rock_obs, bina_worm_obs, bina_wuer_obs,
                  bina_kaub_obs, bina_coch_obs, bina_koel_obs, bina_koel_obs)

bina_sim <- cbind(bina_base_sim, bina_spey_sim, bina_rock_sim, bina_worm_sim, bina_wuer_sim,
                  bina_kaub_sim, bina_coch_sim, bina_koel_sim, bina_koel_sim)

for(i in 1:ncol(bina_obs)){
  
  eca_out <- CC.eca.ts(seriesA = bina_obs[, i], seriesB = bina_sim[, i], delT = 7, tau = 0,
                       sym = TRUE, sigtest = "poisson")[5]
  #precursor coincidence: fraction of A-type events preceded by at least one B-type event
  
  print(eca_out)
  
}


#snow_cover----

#Snow cover duration simulations

nc_flux_file <- paste0(run_dir, "output/mHM_Fluxes_States.nc")
nc_flux <- nc_open(nc_flux_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

sta_date_ind <- which(format(date) == "1954-01-02")
count_date <- length(date)

#Fluxes and states
snow_cube <- ncvar_get(nc_flux, start = c(1, 1, sta_date_ind), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")

sd2sc <- function(val_in, sc_thr = 2){
  
  if(is.na(val_in)){
    val_out <- NA
  }else{
    if(val_in > sc_thr){
      val_out <- 1
    }else{
      val_out <- 0
    }
    
  }
  
}

date_scd <- scf_date

date_scd_ind <- which(date %in% date_scd)

cells_sel <- which(!is.na(c(snow_cube[, , date_scd_ind[1]])))

for(i in 1:length(date_scd_ind)){
  
  print(i)
  
  scd_sim_sing <- sapply(c(snow_cube[, , date_scd_ind[i]][cells_sel]), sd2sc)
  
  if(i == 1){
    scd_sim <-  scd_sim_sing
  }else{
    scd_sim <-  scd_sim + scd_sim_sing
  }
  
}

#Snow cover duration MODIS

#get file names
file_names <- dir(path = scf_dlr_dir, recursive = T)
# unique(nchar(file_names))
# file_names[which(nchar(file_names) %in% c(36, 32, 10))]
file_names <- file_names[which(nchar(file_names) == nchar(file_names[1]))]
scf_file <- raster(paste0(scf_dlr_dir , file_names[1]))

basin_lobi_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/lobith_catch.shp")

basin_lobi_buf <- buffer(basin_lobi_raw, width = 30000)

basin_lobi <- spTransform(basin_lobi_buf, CRS = crs(scf_file, asText = T))

#get dates from file names
f_scf_date <- function(file_path, provider = "DLR"){
  
  #Extract date from file name
  if(provider == "DLR"){
    
    doy <- as.numeric(substr(file_path, nchar(file_path)-6, nchar(file_path)-4))
    yea <- substr(file_path, nchar(file_path)-11, nchar(file_path)-8)
    date <- as.character(as.Date(doy, origin = paste0(yea, "-01-01")))
    
  }
  
  if(provider == "EURAC"){
    
    day <- substr(file_path, nchar(file_path)-12, nchar(file_path)-11)
    mon <- substr(file_path, nchar(file_path)-14, nchar(file_path)-13)
    yea <- substr(file_path, nchar(file_path)-18, nchar(file_path)-15)
    date <- paste0(yea, "-", mon, "-", day)
    
  }
  
  return(date)
  
}

f_scf_date(paste0(scf_dlr_dir, file_names[1]))

scf_dates <- foreach(i = 1:length(file_names), .combine = 'cbind') %dopar% {
  
  f_scf_date(paste0(scf_dlr_dir, file_names[i]))
  
}

scf_date <- as.Date(as.character(scf_dates[1, ]), "%Y-%m-%d")

# #select files covering validation period
# file_names_sel <- file_names[which(scf_date %in% date_vali)]
file_names_sel <- file_names

#sum up days with snow cover for selected basin (with buffer) and time frame

f_scd_extr <- function(file_path, snow_val, basin_in, provider){
  
  #Read file
  scf <- raster(file_path)

  #corp file to basin area (with buffer)
  scf_cro <- raster::crop(scf, extent(basin_in))
  scf_sub <- raster::mask(scf_cro, basin_in)
  
  #get values and set to 0 if not snow, to 1 if snow
  scf_val_NA <- values(scf_sub)
  scf_val_NA[which(scf_val_NA != snow_val)] <- 0
  scf_val_NA[which(scf_val_NA == snow_val)] <- 1
  
  #Extract date from file name
  if(provider == "DLR"){
    
    doy <- as.numeric(substr(file_path, nchar(file_path)-6, nchar(file_path)-4))
    yea <- substr(file_path, nchar(file_path)-11, nchar(file_path)-8)
    date <- as.character(as.Date(doy, origin = paste0(yea, "-01-01")))
    
  }
  
  if(provider == "EURAC"){
    
    day <- substr(file_path, nchar(file_path)-12, nchar(file_path)-11)
    mon <- substr(file_path, nchar(file_path)-14, nchar(file_path)-13)
    yea <- substr(file_path, nchar(file_path)-18, nchar(file_path)-15)
    date <- paste0(yea, "-", mon, "-", day)
    
  }
  
  return(scf_val_NA)
  
}

f_scd_extr(file_path = paste0(scf_dlr_dir, file_names_calc[1]),
           snow_val = 50,
           basin_in = basin_lobi,
           provider = "DLR")

block_size <- 1000
block_stas <- c(1, seq(block_size+1, length(file_names_sel), by = block_size))
block_ends <- c(seq(block_size, length(file_names_sel), by = block_size), length(file_names_sel))

for(b in 1:length(block_stas)){
  
  file_names_calc <- file_names_sel[block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Spatial analysis: Days with snow cover", "Block:", b, "out of", length(block_stas)))
  
  scd_out <- foreach(i = 1:length(file_names_calc), .combine = 'cbind') %dopar% {
    
    f_scd_extr(file_path = paste0(scf_dlr_dir, file_names_calc[i]),
               snow_val = 50,
               basin_in = basin_lobi,
               provider = "DLR")
    
  }
  
  if(b == 1){
    scd_buf_all <- scd_out
  }else{
    scd_buf_all <- cbind(scd_buf_all, scd_out)
  }
}

scd_buf_sum <- apply(scd_buf_all, 1, sum_na)
rm(scd_buf_all) #remove file after calculation as very big
gc() #colltect some garbage

#fill dummy raster with calculated snow cover fraction values
scf_buf_crop <- raster::crop(scf_file, extent(basin_lobi))
scf_buf <- mask(scf_buf_crop, basin_lobi)
scf_buf@data@values <- scd_buf_sum
scf_buf_aggr <- aggregate(scf_buf, fact = 10, fun = mean, na.rm = TRUE)
plot(scf_buf, col = viridis(200, direction = -1))
plot(scf_buf_aggr, col = viridis(200, direction = -1))

#Get values grid points simulated

grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)

scd_dlr <- raster::extract(scf_buf_aggr, grid_points_cube_84[cells_sel])


val2col <- function(val_in, dat_ref, do_log = F, do_bicol = T, col_na = "white"){
  
  if(do_log){
    
    val_in <- log(val_in)
    dat_ref <- log(dat_ref)
    
  }
  
  if(is.na(val_in)){#set NAs to mean to keep script running; later back to NA
    val_in <- mea_na(dat_ref)
    set2NA_1 <- T
  }else{
    set2NA_1 <- F
  }
  
  if(do_bicol){
    
    col_ind <- round((abs(val_in) / max_na(abs(dat_ref))) * 100)
    
    if(val_in < 0){
      my_col  <- colorRampPalette(c("grey80", "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "orangered4", "darkred"))(100)
    }else{
      my_col  <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
    }
    
  }else{
    col_ind <- round((val_in-min_na(dat_ref)) / (max_na(dat_ref)-min_na(dat_ref)) * 200)  
    my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
  }
  
  
  if(is.na(col_ind)){
    set2NA_2 <- T
    col_ind <- 1 #set to one to keep script running; later set to NA color
  }else{
    set2NA_2 = F
  }
  
  if(col_ind == 0){#for minimum and very small values
    
    col_ind <- 1
    
  }
  
  col_out <- my_col[col_ind]
  
  if(length(col_out) < 1){
    
    col_out <- col_na
    
  }
  
  if(set2NA_1 | set2NA_2){
    
    col_out <- col_na
    
  }
  
  return(col_out)
  
}

#Values to colors simulation
scd_sim_ann <- scd_sim / round(length(date_scd) / 365)
cols_spat_sim <- foreach(i = 1:length(scd_sim_ann), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_sim_ann[i],
          dat_ref = scd_sim_ann,
          do_bicol = F)
  
}

#Values to colors difference
scd_dif <- (scd_sim - scd_dlr) / round(length(date_scd) / 365) #Calculate difference Obs. and Sim.
cols_spat_dif <- foreach(i = 1:length(scd_dif), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif[i], 
          dat_ref = scd_dif,
          do_log = F,
          do_bicol = T)
  
}

#Values to colors observations
scd_dlr_ann <- scd_dlr / round(length(date_scd) / 365)
cols_spat_obs <- foreach(i = 1:length(scd_dlr_ann), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dlr_ann[i],
          dat_ref = scd_dlr_ann,
          do_bicol = F)
  
}


pdf(paste0(bas_dir, "res_figs/scd_maps",".pdf"), width = 16, height = 4.2)

#Plot maps
layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4, rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())
# layout.show(n = 7)

par(family = "serif")
cex_pch <- 0.60

#Map Simulations
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_sim, cex = cex_pch)
# plot(basin_base, add =T, lwd = 1.5)
mtext("a) Snow simulations", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(scd_sim_ann)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_sim_ann), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map Difference
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_dif, cex = cex_pch)
# plot(basin_base, add = T, lwd = 1.5)
mtext("b) Difference (a - c)", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(-max_na(abs(scd_dif)), max_na(abs(scd_dif)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map MODIS
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_obs, cex = cex_pch)
# plot(basin_base, add =T, lwd = 1.5)
mtext("c) MODIS snow cover", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(scd_dlr_ann)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dlr_ann), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

dev.off()


#Changes snow cover duration over time

date_scd_1 <- date_simu <- seq(as.Date("1954-01-01", format = "%Y-%m-%d"), 
                               as.Date("1983-12-31", format = "%Y-%m-%d"), by = "day")

date_scd_2 <- date_simu <- seq(as.Date("1984-01-01", format = "%Y-%m-%d"), 
                               as.Date("2013-12-31", format = "%Y-%m-%d"), by = "day")

date_scd_ind_1 <- which(date %in% date_scd_1)
date_scd_ind_2 <- which(date %in% date_scd_2)

cells_sel <- which(!is.na(c(snow_cube[, , date_scd_ind[1]])))

for(i in 1:length(date_scd_ind_1)){
  
  print(i)
  
  scd_sim_sing <- sapply(c(snow_cube[, , date_scd_ind_1[i]][cells_sel]), sd2sc)
  
  if(i == 1){
    scd_sim_1 <-  scd_sim_sing
  }else{
    scd_sim_1 <-  scd_sim_1 + scd_sim_sing
  }
  
}

for(i in 1:length(date_scd_ind_2)){
  
  print(i)
  
  scd_sim_sing <- sapply(c(snow_cube[, , date_scd_ind_2[i]][cells_sel]), sd2sc)
  
  if(i == 1){
    scd_sim_2 <-  scd_sim_sing
  }else{
    scd_sim_2 <-  scd_sim_2 + scd_sim_sing
  }
  
}

#Values to colors simulation part 1
scd_sim_ann_1 <- scd_sim_1 / round(length(date_scd_1) / 365)
cols_spat_sim_1 <- foreach(i = 1:length(scd_sim_ann_1), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_sim_ann_1[i],
          dat_ref = scd_sim_ann_1,
          do_bicol = F)
  
}

#Values to colors difference
scd_dif <- (scd_sim_1 - scd_sim_2) / round(length(date_scd_1) / 365) #Calculate difference Obs. and Sim.
cols_spat_dif <- foreach(i = 1:length(scd_dif), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif[i], 
          dat_ref = scd_dif,
          do_log = F,
          do_bicol = T)
  
}

#Values to colors observations
scd_sim_ann_2 <- scd_sim_2 / round(length(date_scd_2) / 365)
cols_spat_obs <- foreach(i = 1:length(scd_sim_ann_2), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_sim_ann_2[i],
          dat_ref = scd_sim_ann_2,
          do_bicol = F)
  
}


pdf(paste0(bas_dir, "res_figs/scd_maps_time",".pdf"), width = 16, height = 4.2)

#Plot maps
layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4, rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())
# layout.show(n = 7)

par(family = "serif")
cex_pch <- 0.60

#Map Simulations
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_sim, cex = cex_pch)
# plot(basin_base, add =T, lwd = 1.5)
mtext("a) Simulations 1954-1983", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(scd_sim_ann_1)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_sim_ann_1), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map Difference
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_dif, cex = cex_pch)
# plot(basin_base, add = T, lwd = 1.5)
mtext("b) Difference (a - c)", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(-max_na(abs(scd_dif)), max_na(abs(scd_dif)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map MODIS
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], pch = 15, col = cols_spat_obs, cex = cex_pch)
# plot(basin_base, add =T, lwd = 1.5)
mtext("c) Simulations 1984-2013", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(scd_sim_ann_2)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_sim_ann_2), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

dev.off()
