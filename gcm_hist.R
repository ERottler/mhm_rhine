###

#mHM simulations Part II: mHM simulations using historical GCM data
#Erwin Rottler, University of Potsdam, Summer 2020

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales, ncdf4,
               ncdf4.helpers, sp, raster, viridis, meltimr, POT, readr, hydroGOF,
               CoinCalc, seas)

#set directories
bas_dir <- "U:/rhine_fut/R/"
run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

#load functions
source(paste0(bas_dir, "mhm_rhine/functs.R"))

#evaluation time frame
sta_yea <- 1951
end_yea <- 2000

date_simu <- seq(as.Date("1951-01-01", format = "%Y-%m-%d"), 
                 as.Date("2000-12-31", format = "%Y-%m-%d"), by = "day")

stopCluster(my_clust)

n_cores <- 5 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#get_obs_disc----

#get meta data and measured time series for selected gauges (GRDC)

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

#get_sim_disc----

nc_disc_file_obs <- paste0(run_dir, "output/EOBS/output/mRM_Fluxes_States.nc")
nc_disc_file_cm1 <- paste0(run_dir, "output/GFDL-ESM2M/historical/output/mRM_Fluxes_States.nc")
nc_disc_file_cm2 <- paste0(run_dir, "output/HadGEM2-ES/historical/output/mRM_Fluxes_States.nc")
nc_disc_file_cm3 <- paste0(run_dir, "output/IPSL-CM5A-LR/historical/output/mRM_Fluxes_States.nc")
nc_disc_file_cm4 <- paste0(run_dir, "output/MIROC-ESM-CHEM/historical/output/mRM_Fluxes_States.nc")
nc_disc_file_cm5 <- paste0(run_dir, "output/NorESM1-M/historical/output/mRM_Fluxes_States.nc")

nc_disc_obs <- ncdf4::nc_open(nc_disc_file_obs)
nc_disc_cm1 <- ncdf4::nc_open(nc_disc_file_cm1)
nc_disc_cm2 <- ncdf4::nc_open(nc_disc_file_cm2)
nc_disc_cm3 <- ncdf4::nc_open(nc_disc_file_cm3)
nc_disc_cm4 <- ncdf4::nc_open(nc_disc_file_cm4)
nc_disc_cm5 <- ncdf4::nc_open(nc_disc_file_cm5)

#Get simulated runoff for selected gauges
lon <- ncdf4::ncvar_get(nc_disc_obs, varid = "lon")
lat <- ncdf4::ncvar_get(nc_disc_obs, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_disc_obs, time.dim.name = "time")))

disc_cube <- ncvar_get(nc_disc_obs, start = c(1, 1, 1),
                       count = c(nrow(lon), ncol(lon), 1), varid = "Qrouted")

gaug_spa <- SpatialPoints(grdc_meta[, c(2, 3)], proj4string = crswgs84)
grid_spa <- SpatialPoints(cbind(c(lat), c(lon)), proj4string = crswgs84)

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

rows_sel_gaugs[3] <- rows_sel_gaugs[3]+1 #gauge cochem one row lower


disc_from_nc <- function(nc_disc, date_sel = date_simu, varID = "Qrouted"){
  
  date <- as.Date(as.character(nc.get.time.series(nc_disc, time.dim.name = "time")))
  count_date <- length(date)
  
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
  
  disc_sim_full <- cbind(simu_lobi, simu_koel, simu_coch, simu_kaub, simu_wuer, simu_worm, simu_rock, 
                    simu_spey, simu_base, simu_unte, simu_reki)
  
  disc_sim <- disc_sim_full[which(date %in% date_sel), ]
  
  return(disc_sim)
  
}

disc_obs <- as.data.frame(disc_from_nc(nc_disc_obs))
disc_cm1 <- as.data.frame(disc_from_nc(nc_disc_cm1))
disc_cm2 <- as.data.frame(disc_from_nc(nc_disc_cm2))
disc_cm3 <- as.data.frame(disc_from_nc(nc_disc_cm3))
disc_cm4 <- as.data.frame(disc_from_nc(nc_disc_cm4))
disc_cm5 <- as.data.frame(disc_from_nc(nc_disc_cm5))

#high_flow_indicators----

f_ann_max <- function(data_in, date = date_simu, start_y = 1951, end_y = 2000, break_day = 274){

  disc_day <- meltimr::ord_day(data_in = data_in,
                               date = date,
                               start_y = start_y,
                               end_y = end_y,
                               break_day = break_day)
  
  disc_max <- apply(disc_day, 1, max_na)
  
  return(disc_max)
  
}
f_ann_q90 <- function(data_in, date = date_simu, start_y = 1951, end_y = 2000, break_day = 274){
  
  disc_day <- meltimr::ord_day(data_in = data_in,
                               date = date,
                               start_y = start_y,
                               end_y = end_y,
                               break_day = break_day)
  
  q90 <- function(data_in){
    
    stats::quantile(x = data_in, probs = c(0.9), na.rm = TRUE, type = 8)
    
  }
  
  disc_q90 <- apply(disc_day, 1, q90)
  
  return(disc_q90)
  
}

#Observations/Measurements from GRDC
disc_mea_max <- cbind(f_ann_max(disc_lobi$value), f_ann_max(disc_koel$value), f_ann_max(disc_coch$value),
                      f_ann_max(disc_kaub$value), f_ann_max(disc_wuer$value), f_ann_max(disc_worm$value),
                      f_ann_max(disc_rock$value), f_ann_max(disc_spey$value), f_ann_max(disc_base$value), 
                      f_ann_max(disc_unte$value), f_ann_max(disc_reki$value))

disc_mea_q90 <- cbind(f_ann_q90(disc_lobi$value), f_ann_q90(disc_koel$value), f_ann_q90(disc_coch$value),
                      f_ann_q90(disc_kaub$value), f_ann_q90(disc_wuer$value), f_ann_q90(disc_worm$value),
                      f_ann_q90(disc_rock$value), f_ann_q90(disc_spey$value), f_ann_q90(disc_base$value), 
                      f_ann_q90(disc_unte$value), f_ann_q90(disc_reki$value))

#Simulations based on observations
disc_obs_max <- cbind(f_ann_max(disc_obs$simu_lobi), f_ann_max(disc_obs$simu_koel), f_ann_max(disc_obs$simu_coch),
                      f_ann_max(disc_obs$simu_kaub), f_ann_max(disc_obs$simu_wuer), f_ann_max(disc_obs$simu_worm),
                      f_ann_max(disc_obs$simu_rock), f_ann_max(disc_obs$simu_spey), f_ann_max(disc_obs$simu_base), 
                      f_ann_max(disc_obs$simu_unte), f_ann_max(disc_obs$simu_reki))

disc_obs_q90 <- cbind(f_ann_q90(disc_obs$simu_lobi), f_ann_q90(disc_obs$simu_koel), f_ann_q90(disc_obs$simu_coch),
                      f_ann_q90(disc_obs$simu_kaub), f_ann_q90(disc_obs$simu_wuer), f_ann_q90(disc_obs$simu_worm),
                      f_ann_q90(disc_obs$simu_rock), f_ann_q90(disc_obs$simu_spey), f_ann_q90(disc_obs$simu_base), 
                      f_ann_q90(disc_obs$simu_unte), f_ann_q90(disc_obs$simu_reki))

#Simulations based on cm1
disc_cm1_max <- cbind(f_ann_max(disc_cm1$simu_lobi), f_ann_max(disc_cm1$simu_koel), f_ann_max(disc_cm1$simu_coch),
                      f_ann_max(disc_cm1$simu_kaub), f_ann_max(disc_cm1$simu_wuer), f_ann_max(disc_cm1$simu_worm),
                      f_ann_max(disc_cm1$simu_rock), f_ann_max(disc_cm1$simu_spey), f_ann_max(disc_cm1$simu_base), 
                      f_ann_max(disc_cm1$simu_unte), f_ann_max(disc_cm1$simu_reki))

disc_cm1_q90 <- cbind(f_ann_q90(disc_cm1$simu_lobi), f_ann_q90(disc_cm1$simu_koel), f_ann_q90(disc_cm1$simu_coch),
                      f_ann_q90(disc_cm1$simu_kaub), f_ann_q90(disc_cm1$simu_wuer), f_ann_q90(disc_cm1$simu_worm),
                      f_ann_q90(disc_cm1$simu_rock), f_ann_q90(disc_cm1$simu_spey), f_ann_q90(disc_cm1$simu_base), 
                      f_ann_q90(disc_cm1$simu_unte), f_ann_q90(disc_cm1$simu_reki))

#Simulations based on cm2
disc_cm2_max <- cbind(f_ann_max(disc_cm2$simu_lobi), f_ann_max(disc_cm2$simu_koel), f_ann_max(disc_cm2$simu_coch),
                      f_ann_max(disc_cm2$simu_kaub), f_ann_max(disc_cm2$simu_wuer), f_ann_max(disc_cm2$simu_worm),
                      f_ann_max(disc_cm2$simu_rock), f_ann_max(disc_cm2$simu_spey), f_ann_max(disc_cm2$simu_base), 
                      f_ann_max(disc_cm2$simu_unte), f_ann_max(disc_cm2$simu_reki))

disc_cm2_q90 <- cbind(f_ann_q90(disc_cm2$simu_lobi), f_ann_q90(disc_cm2$simu_koel), f_ann_q90(disc_cm2$simu_coch),
                      f_ann_q90(disc_cm2$simu_kaub), f_ann_q90(disc_cm2$simu_wuer), f_ann_q90(disc_cm2$simu_worm),
                      f_ann_q90(disc_cm2$simu_rock), f_ann_q90(disc_cm2$simu_spey), f_ann_q90(disc_cm2$simu_base), 
                      f_ann_q90(disc_cm2$simu_unte), f_ann_q90(disc_cm2$simu_reki))

#Simulations based on cm3
disc_cm3_max <- cbind(f_ann_max(disc_cm3$simu_lobi), f_ann_max(disc_cm3$simu_koel), f_ann_max(disc_cm3$simu_coch),
                      f_ann_max(disc_cm3$simu_kaub), f_ann_max(disc_cm3$simu_wuer), f_ann_max(disc_cm3$simu_worm),
                      f_ann_max(disc_cm3$simu_rock), f_ann_max(disc_cm3$simu_spey), f_ann_max(disc_cm3$simu_base), 
                      f_ann_max(disc_cm3$simu_unte), f_ann_max(disc_cm3$simu_reki))

disc_cm3_q90 <- cbind(f_ann_q90(disc_cm3$simu_lobi), f_ann_q90(disc_cm3$simu_koel), f_ann_q90(disc_cm3$simu_coch),
                      f_ann_q90(disc_cm3$simu_kaub), f_ann_q90(disc_cm3$simu_wuer), f_ann_q90(disc_cm3$simu_worm),
                      f_ann_q90(disc_cm3$simu_rock), f_ann_q90(disc_cm3$simu_spey), f_ann_q90(disc_cm3$simu_base), 
                      f_ann_q90(disc_cm3$simu_unte), f_ann_q90(disc_cm3$simu_reki))

#Simulations based on cm4
disc_cm4_max <- cbind(f_ann_max(disc_cm4$simu_lobi), f_ann_max(disc_cm4$simu_koel), f_ann_max(disc_cm4$simu_coch),
                      f_ann_max(disc_cm4$simu_kaub), f_ann_max(disc_cm4$simu_wuer), f_ann_max(disc_cm4$simu_worm),
                      f_ann_max(disc_cm4$simu_rock), f_ann_max(disc_cm4$simu_spey), f_ann_max(disc_cm4$simu_base), 
                      f_ann_max(disc_cm4$simu_unte), f_ann_max(disc_cm4$simu_reki))

disc_cm4_q90 <- cbind(f_ann_q90(disc_cm4$simu_lobi), f_ann_q90(disc_cm4$simu_koel), f_ann_q90(disc_cm4$simu_coch),
                      f_ann_q90(disc_cm4$simu_kaub), f_ann_q90(disc_cm4$simu_wuer), f_ann_q90(disc_cm4$simu_worm),
                      f_ann_q90(disc_cm4$simu_rock), f_ann_q90(disc_cm4$simu_spey), f_ann_q90(disc_cm4$simu_base), 
                      f_ann_q90(disc_cm4$simu_unte), f_ann_q90(disc_cm4$simu_reki))

#Simulations based on cm5
disc_cm5_max <- cbind(f_ann_max(disc_cm5$simu_lobi), f_ann_max(disc_cm5$simu_koel), f_ann_max(disc_cm5$simu_coch),
                      f_ann_max(disc_cm5$simu_kaub), f_ann_max(disc_cm5$simu_wuer), f_ann_max(disc_cm5$simu_worm),
                      f_ann_max(disc_cm5$simu_rock), f_ann_max(disc_cm5$simu_spey), f_ann_max(disc_cm5$simu_base), 
                      f_ann_max(disc_cm5$simu_unte), f_ann_max(disc_cm5$simu_reki))

disc_cm5_q90 <- cbind(f_ann_q90(disc_cm5$simu_lobi), f_ann_q90(disc_cm5$simu_koel), f_ann_q90(disc_cm5$simu_coch),
                      f_ann_q90(disc_cm5$simu_kaub), f_ann_q90(disc_cm5$simu_wuer), f_ann_q90(disc_cm5$simu_worm),
                      f_ann_q90(disc_cm5$simu_rock), f_ann_q90(disc_cm5$simu_spey), f_ann_q90(disc_cm5$simu_base), 
                      f_ann_q90(disc_cm5$simu_unte), f_ann_q90(disc_cm5$simu_reki))

lims_max <- c(1, max_na(c(disc_mea_max, disc_obs_max, disc_cm1_max, disc_cm1_max, disc_cm2_max, 
                          disc_cm3_max, disc_cm4_max, disc_cm5_max)))
lims_q90 <- c(1, max_na(c(disc_mea_q90, disc_obs_q90, disc_cm1_q90, disc_cm1_q90, disc_cm2_q90, 
                          disc_cm3_q90, disc_cm4_q90, disc_cm5_q90)))

#plot_eval----

pdf(paste0(bas_dir,"res_figs/eval_hist.pdf"), width = 16, height = 6)

par(family = "serif")
par(mar = c(3.5, 3.5, 2.5, 0.5))
par(mfrow = c(2, 4))
cex_points <- 1.4

#Meased vs. EOBS Q90 all stations
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_obs_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_obs_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("Q90 Measured vs. EOBS", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("All gauges validation"), pch = 19, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Meased vs. EOBS MAX all stations
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_obs_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_obs_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("MAX Measured vs. EOBS", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("All gauges validation"), pch = 19, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Meased vs. EOBS Q90 selected stations
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_obs_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_obs_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("Q90 Measured vs. EOBS", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("Cologne", "Basel", "Cochem"), pch = 19, 
       col = c("black", "darkred", "steelblue4"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Meased vs. EOBS Max selected stations
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_obs_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_obs_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("Q90 Measured vs. EOBS", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("Cologne", "Basel", "Cochem"), pch = 19, 
       col = c("black", "darkred", "steelblue4"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Measured vs. GCMs Q90
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_mea_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_cm1_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm2_q90[, i], pch = 22, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm3_q90[, i], pch = 23, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm4_q90[, i], pch = 24, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm5_q90[, i], pch = 25, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)

}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("Q90 Measured vs. GCMs", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("All gauges validation"), pch = 19, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
legend("bottomright", c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM", "NorESM1-M"), pch = 21:25, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Measured vs. GCMs MAX
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_mea_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_cm1_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm2_max[, i], pch = 22, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm3_max[, i], pch = 23, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm4_max[, i], pch = 24, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm5_max[, i], pch = 25, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("MAX Measured vs. GCMs", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("All gauges validation"), pch = 19, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
legend("bottomright", c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM", "NorESM1-M"), pch = 21:25, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Measured vs. GCMs selected Q90
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_mea_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_cm1_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm2_q90[, i], pch = 22, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm3_q90[, i], pch = 23, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm4_q90[, i], pch = 24, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_q90[, i], disc_cm5_q90[, i], pch = 25, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("Q90 Measured vs. GCMs", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("Cologne", "Basel", "Cochem"), pch = 19, 
       col = c("black", "darkred", "steelblue4"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
legend("bottomright", c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM", "NorESM1-M"), pch = 21:25, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Measured vs. GCMs selected MAX
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_mea_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_cm1_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm2_max[, i], pch = 22, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm3_max[, i], pch = 23, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm4_max[, i], pch = 24, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  points(disc_mea_max[, i], disc_cm5_max[, i], pch = 25, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
abline(a = 0, b = 1)
axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
axis(2, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.3)
mtext("MAX Measured vs. GCMs", side = 3, line = 0.5, cex = 1.3)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 1, line = 1.8, cex = 1.1)
mtext(expression(paste("Discharge [m"^"3", "s"^"-1","]")), side = 2, line = 1.4, cex = 1.1)
legend("topleft", c("Cologne", "Basel", "Cochem"), pch = 19, 
       col = c("black", "darkred", "steelblue4"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
legend("bottomright", c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM", "NorESM1-M"), pch = 21:25, 
       col = c("black"), cex = 1.0,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

dev.off()









