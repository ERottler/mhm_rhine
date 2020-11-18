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
# unte_file <- paste0(grdc_dir, "6935300_Q_Day.Cmd.txt")
# reki_file <- paste0(grdc_dir, "6935054_Q_Day.Cmd.txt")

file_paths <- c(lobi_file, koel_file, coch_file, kaub_file, wuer_file, worm_file,
                rock_file, spey_file, base_file)

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

disc_lobi <- disc_lobi_full[which(disc_lobi_full$date %in% date_simu), ]
disc_koel <- disc_koel_full[which(disc_koel_full$date %in% date_simu), ]
disc_coch <- disc_coch_full[which(disc_coch_full$date %in% date_simu), ]
disc_kaub <- disc_kaub_full[which(disc_kaub_full$date %in% date_simu), ]
disc_wuer <- disc_wuer_full[which(disc_wuer_full$date %in% date_simu), ]
disc_worm <- disc_worm_full[which(disc_worm_full$date %in% date_simu), ]
disc_rock <- disc_rock_full[which(disc_rock_full$date %in% date_simu), ]
disc_spey <- disc_spey_full[which(disc_spey_full$date %in% date_simu), ]
disc_base <- disc_base_full[which(disc_base_full$date %in% date_simu), ]

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
  
  disc_sim_full <- cbind(simu_lobi, simu_koel, simu_coch, simu_kaub, simu_wuer, simu_worm, simu_rock, 
                         simu_spey, simu_base)
  
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
f_ann_doy <- function(data_in, date = date_simu, start_y = 1951, end_y = 2000, break_day = 274){
  
  disc_day <- meltimr::ord_day(data_in = data_in,
                               date = date,
                               start_y = start_y,
                               end_y = end_y,
                               break_day = break_day)
  
  doy_max <- function(data_in){
    
    min_na(which(data_in == max_na(data_in)))
    
  }
  
  disc_doy <- apply(disc_day, 1, doy_max)
  
  return(disc_doy)
  
}

#Observations/Measurements from GRDC
disc_mea_max <- cbind(f_ann_max(disc_lobi$value), f_ann_max(disc_koel$value), f_ann_max(disc_coch$value),
                      f_ann_max(disc_kaub$value), f_ann_max(disc_wuer$value), f_ann_max(disc_worm$value),
                      f_ann_max(disc_rock$value), f_ann_max(disc_spey$value), f_ann_max(disc_base$value))

disc_mea_q90 <- cbind(f_ann_q90(disc_lobi$value), f_ann_q90(disc_koel$value), f_ann_q90(disc_coch$value),
                      f_ann_q90(disc_kaub$value), f_ann_q90(disc_wuer$value), f_ann_q90(disc_worm$value),
                      f_ann_q90(disc_rock$value), f_ann_q90(disc_spey$value), f_ann_q90(disc_base$value))

disc_mea_doy <- cbind(f_ann_doy(disc_lobi$value), f_ann_doy(disc_koel$value), f_ann_doy(disc_coch$value),
                      f_ann_doy(disc_kaub$value), f_ann_doy(disc_wuer$value), f_ann_doy(disc_worm$value),
                      f_ann_doy(disc_rock$value), f_ann_doy(disc_spey$value), f_ann_doy(disc_base$value))

#Simulations based on observations
disc_obs_max <- cbind(f_ann_max(disc_obs$simu_lobi), f_ann_max(disc_obs$simu_koel), f_ann_max(disc_obs$simu_coch),
                      f_ann_max(disc_obs$simu_kaub), f_ann_max(disc_obs$simu_wuer), f_ann_max(disc_obs$simu_worm),
                      f_ann_max(disc_obs$simu_rock), f_ann_max(disc_obs$simu_spey), f_ann_max(disc_obs$simu_base))

disc_obs_q90 <- cbind(f_ann_q90(disc_obs$simu_lobi), f_ann_q90(disc_obs$simu_koel), f_ann_q90(disc_obs$simu_coch),
                      f_ann_q90(disc_obs$simu_kaub), f_ann_q90(disc_obs$simu_wuer), f_ann_q90(disc_obs$simu_worm),
                      f_ann_q90(disc_obs$simu_rock), f_ann_q90(disc_obs$simu_spey), f_ann_q90(disc_obs$simu_base))

disc_obs_doy <- cbind(f_ann_doy(disc_obs$simu_lobi), f_ann_doy(disc_obs$simu_koel), f_ann_doy(disc_obs$simu_coch),
                      f_ann_doy(disc_obs$simu_kaub), f_ann_doy(disc_obs$simu_wuer), f_ann_doy(disc_obs$simu_worm),
                      f_ann_doy(disc_obs$simu_rock), f_ann_doy(disc_obs$simu_spey), f_ann_doy(disc_obs$simu_base))

#Simulations based on cm1
disc_cm1_max <- cbind(f_ann_max(disc_cm1$simu_lobi), f_ann_max(disc_cm1$simu_koel), f_ann_max(disc_cm1$simu_coch),
                      f_ann_max(disc_cm1$simu_kaub), f_ann_max(disc_cm1$simu_wuer), f_ann_max(disc_cm1$simu_worm),
                      f_ann_max(disc_cm1$simu_rock), f_ann_max(disc_cm1$simu_spey), f_ann_max(disc_cm1$simu_base))

disc_cm1_q90 <- cbind(f_ann_q90(disc_cm1$simu_lobi), f_ann_q90(disc_cm1$simu_koel), f_ann_q90(disc_cm1$simu_coch),
                      f_ann_q90(disc_cm1$simu_kaub), f_ann_q90(disc_cm1$simu_wuer), f_ann_q90(disc_cm1$simu_worm),
                      f_ann_q90(disc_cm1$simu_rock), f_ann_q90(disc_cm1$simu_spey), f_ann_q90(disc_cm1$simu_base))

disc_cm1_doy <- cbind(f_ann_doy(disc_cm1$simu_lobi), f_ann_doy(disc_cm1$simu_koel), f_ann_doy(disc_cm1$simu_coch),
                      f_ann_doy(disc_cm1$simu_kaub), f_ann_doy(disc_cm1$simu_wuer), f_ann_doy(disc_cm1$simu_worm),
                      f_ann_doy(disc_cm1$simu_rock), f_ann_doy(disc_cm1$simu_spey), f_ann_doy(disc_cm1$simu_base))

#Simulations based on cm2
disc_cm2_max <- cbind(f_ann_max(disc_cm2$simu_lobi), f_ann_max(disc_cm2$simu_koel), f_ann_max(disc_cm2$simu_coch),
                      f_ann_max(disc_cm2$simu_kaub), f_ann_max(disc_cm2$simu_wuer), f_ann_max(disc_cm2$simu_worm),
                      f_ann_max(disc_cm2$simu_rock), f_ann_max(disc_cm2$simu_spey), f_ann_max(disc_cm2$simu_base))

disc_cm2_q90 <- cbind(f_ann_q90(disc_cm2$simu_lobi), f_ann_q90(disc_cm2$simu_koel), f_ann_q90(disc_cm2$simu_coch),
                      f_ann_q90(disc_cm2$simu_kaub), f_ann_q90(disc_cm2$simu_wuer), f_ann_q90(disc_cm2$simu_worm),
                      f_ann_q90(disc_cm2$simu_rock), f_ann_q90(disc_cm2$simu_spey), f_ann_q90(disc_cm2$simu_base))

disc_cm2_doy <- cbind(f_ann_doy(disc_cm2$simu_lobi), f_ann_doy(disc_cm2$simu_koel), f_ann_doy(disc_cm2$simu_coch),
                      f_ann_doy(disc_cm2$simu_kaub), f_ann_doy(disc_cm2$simu_wuer), f_ann_doy(disc_cm2$simu_worm),
                      f_ann_doy(disc_cm2$simu_rock), f_ann_doy(disc_cm2$simu_spey), f_ann_doy(disc_cm2$simu_base))

#Simulations based on cm3
disc_cm3_max <- cbind(f_ann_max(disc_cm3$simu_lobi), f_ann_max(disc_cm3$simu_koel), f_ann_max(disc_cm3$simu_coch),
                      f_ann_max(disc_cm3$simu_kaub), f_ann_max(disc_cm3$simu_wuer), f_ann_max(disc_cm3$simu_worm),
                      f_ann_max(disc_cm3$simu_rock), f_ann_max(disc_cm3$simu_spey), f_ann_max(disc_cm3$simu_base))

disc_cm3_q90 <- cbind(f_ann_q90(disc_cm3$simu_lobi), f_ann_q90(disc_cm3$simu_koel), f_ann_q90(disc_cm3$simu_coch),
                      f_ann_q90(disc_cm3$simu_kaub), f_ann_q90(disc_cm3$simu_wuer), f_ann_q90(disc_cm3$simu_worm),
                      f_ann_q90(disc_cm3$simu_rock), f_ann_q90(disc_cm3$simu_spey), f_ann_q90(disc_cm3$simu_base))

disc_cm3_doy <- cbind(f_ann_doy(disc_cm3$simu_lobi), f_ann_doy(disc_cm3$simu_koel), f_ann_doy(disc_cm3$simu_coch),
                      f_ann_doy(disc_cm3$simu_kaub), f_ann_doy(disc_cm3$simu_wuer), f_ann_doy(disc_cm3$simu_worm),
                      f_ann_doy(disc_cm3$simu_rock), f_ann_doy(disc_cm3$simu_spey), f_ann_doy(disc_cm3$simu_base))

#Simulations based on cm4
disc_cm4_max <- cbind(f_ann_max(disc_cm4$simu_lobi), f_ann_max(disc_cm4$simu_koel), f_ann_max(disc_cm4$simu_coch),
                      f_ann_max(disc_cm4$simu_kaub), f_ann_max(disc_cm4$simu_wuer), f_ann_max(disc_cm4$simu_worm),
                      f_ann_max(disc_cm4$simu_rock), f_ann_max(disc_cm4$simu_spey), f_ann_max(disc_cm4$simu_base))

disc_cm4_q90 <- cbind(f_ann_q90(disc_cm4$simu_lobi), f_ann_q90(disc_cm4$simu_koel), f_ann_q90(disc_cm4$simu_coch),
                      f_ann_q90(disc_cm4$simu_kaub), f_ann_q90(disc_cm4$simu_wuer), f_ann_q90(disc_cm4$simu_worm),
                      f_ann_q90(disc_cm4$simu_rock), f_ann_q90(disc_cm4$simu_spey), f_ann_q90(disc_cm4$simu_base))

disc_cm4_doy <- cbind(f_ann_doy(disc_cm4$simu_lobi), f_ann_doy(disc_cm4$simu_koel), f_ann_doy(disc_cm4$simu_coch),
                      f_ann_doy(disc_cm4$simu_kaub), f_ann_doy(disc_cm4$simu_wuer), f_ann_doy(disc_cm4$simu_worm),
                      f_ann_doy(disc_cm4$simu_rock), f_ann_doy(disc_cm4$simu_spey), f_ann_doy(disc_cm4$simu_base))

#Simulations based on cm5
disc_cm5_max <- cbind(f_ann_max(disc_cm5$simu_lobi), f_ann_max(disc_cm5$simu_koel), f_ann_max(disc_cm5$simu_coch),
                      f_ann_max(disc_cm5$simu_kaub), f_ann_max(disc_cm5$simu_wuer), f_ann_max(disc_cm5$simu_worm),
                      f_ann_max(disc_cm5$simu_rock), f_ann_max(disc_cm5$simu_spey), f_ann_max(disc_cm5$simu_base))

disc_cm5_q90 <- cbind(f_ann_q90(disc_cm5$simu_lobi), f_ann_q90(disc_cm5$simu_koel), f_ann_q90(disc_cm5$simu_coch),
                      f_ann_q90(disc_cm5$simu_kaub), f_ann_q90(disc_cm5$simu_wuer), f_ann_q90(disc_cm5$simu_worm),
                      f_ann_q90(disc_cm5$simu_rock), f_ann_q90(disc_cm5$simu_spey), f_ann_q90(disc_cm5$simu_base))

disc_cm5_doy <- cbind(f_ann_doy(disc_cm5$simu_lobi), f_ann_doy(disc_cm5$simu_koel), f_ann_doy(disc_cm5$simu_coch),
                      f_ann_doy(disc_cm5$simu_kaub), f_ann_doy(disc_cm5$simu_wuer), f_ann_doy(disc_cm5$simu_worm),
                      f_ann_doy(disc_cm5$simu_rock), f_ann_doy(disc_cm5$simu_spey), f_ann_doy(disc_cm5$simu_base))

lims_max <- c(10, max_na(c(disc_mea_max, disc_obs_max, disc_cm1_max, disc_cm1_max, disc_cm2_max, 
                          disc_cm3_max, disc_cm4_max, disc_cm5_max)))
lims_q90 <- c(10, max_na(c(disc_mea_q90, disc_obs_q90, disc_cm1_q90, disc_cm1_q90, disc_cm2_q90, 
                          disc_cm3_q90, disc_cm4_q90, disc_cm5_q90)))

#plot_magnitudes----

pdf(paste0(bas_dir,"res_figs/eval_hist_mag.pdf"), width = 16, height = 8)

layout(matrix(c(9, 9, 9, 9,
                1, 2, 3, 4,
                10, 10, 10, 10,
                5, 6, 7, 8),
              4, 4, byrow = T), widths=c(), heights=c(0.15, 1, 0.15, 1))

par(family = "serif")
par(mar = c(4.2, 4.2, 2.5, 1.0))
cex_points <- 1.9
cex_main <- 1.7
main_line <- 0.4
cex_axis <- 1.8
cex_label <- 1.7
x_axs_dist <- 0.5
axs_tic <- 0.012
x_lab_line <- 2.8
y_lab_line <- 1.9
  
#All stations: Q90 Measured  vs. EOBS
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_obs_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_obs_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("a) Q90 Observed vs. EOBS", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[EOBS], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
legend("topleft", c("All gauges"), pch = 19, 
       col = c("black"), cex = 1.8,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#All stations: MAX Measured vs. EOBS
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_obs_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_obs_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("b) MAX Observed vs. EOBS", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[EOBS], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
box()

#All stations: Q90 Measured vs. GCMs
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_mea_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
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
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("c) Q90 Observed vs. GCMs", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[GCM], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
box()

#All stations: MAX Measured vs. GCMs
alpha_sel <- c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
               0.4, 0.4, 0.4, 0.4, 0.4)
cols_sel <- scales::alpha(c("black", "black", "black", "black", "black", "black",
                            "black", "black", "black", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_mea_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
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
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("d) MAX Observed vs. GCMs", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[GCM], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
box()

#Selected stations: Q90 Measured vs. EOBS
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_obs_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
for(i in 1:ncol(disc_mea_q90)){
  
  points(disc_mea_q90[, i], disc_obs_q90[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("e) Q90 Observed vs. EOBS", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[EOBS], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
legend("topleft", c("Cologne", "Basel", "Cochem"), pch = 19, 
       col = c("black", "darkred", "steelblue4"), cex = 1.8,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Selected stations: MAX Measured vs. EOBS
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_obs_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
for(i in 1:ncol(disc_mea_max)){
  
  points(disc_mea_max[, i], disc_obs_max[, i], pch = 21, bg = cols_sel[i], 
         col = alpha("black", alpha = 0), cex = cex_points)
  
}
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("f) MAX Observed vs. EOBS", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[EOBS], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
box()

#Selected stations: Q90 Measured vs. GCMs
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_q90[, 1], disc_mea_q90[, 1], type = "n", log = "xy", ylim = lims_q90, xlim = lims_q90, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
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
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("g) Q90 Observed vs. GCMs", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[GCM], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
box()

#Selectes stations: MAX Measured vs. GCMs
alpha_sel <- c(0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.4, 0.0, 0.0)
cols_sel <- scales::alpha(c("black", "black", "steelblue4", "black", "steelblue4", "black",
                            "black", "black", "darkred", "black", "black"), 
                          alpha = alpha_sel)

plot(disc_mea_max[, 1], disc_mea_max[, 1], type = "n", log = "xy", ylim = lims_max, xlim = lims_max, 
     ylab = "", xlab = "", axes = FALSE)
abline(a = 0, b = 1)
abline(h = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
abline(v = c(50, 500, 5000), lty = "dashed", col = "grey75", lwd = 0.7)
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
axis(1, mgp=c(3, x_axs_dist, 0), tck = -axs_tic, cex.axis = cex_axis)
axis(2, mgp=c(3, 0.19, 0), tck = -axs_tic, cex.axis = cex_axis)
mtext("h) MAX Observed vs. GCMs", side = 3, line = main_line, cex = cex_main, adj = 0)
mtext(expression(paste("Q"[observed], " [m"^"3", "s"^"-1","]")), side = 1, 
      line = x_lab_line, cex = cex_label)
mtext(expression(paste("Q"[GCM], " [m"^"3", "s"^"-1","]")), side = 2, 
      line = y_lab_line, cex = cex_label)
legend("bottomright", c("GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM", "NorESM1-M"), pch = 21:25, 
       col = c("black"), cex = 1.4,
       box.lwd = 0.0, box.col = "black", bg = "white")
box()

#Main header
par(mar = c(0, 0, 0, 0))
plot(1:10, 1:10, type = "n", xlab = "", ylab = "", axes = F)
mtext("All validation gauges", side = 3, line = -2.8, cex = 2.2)

plot(1:10, 1:10, type = "n", xlab = "", ylab = "", axes = F)
mtext("Selected gauges", side = 3, line = -2.8, cex = 2.2)

dev.off()











#plot_timing----

pdf(paste0(bas_dir,"res_figs/eval_hist_doy_raw.pdf"), width = 16, height = 8)

layout(matrix(c(10, 1, 2, 3,
                10, 4, 5, 6,
                10, 7, 8, 9),
              3, 4, byrow = T), widths=c(0.08, 1, 1, 1), heights=c())

par(family = "serif")
par(mar = c(1.0, 4.2, 3.0, 1.0))

plot_doy <- function(doy_mea, doy_obs, doy_cm1, doy_cm2, doy_cm3, doy_cm4, doy_cm5, 
                     calc_ylims = F, ylims = c(0, 365),
                     y_lab = "", do_legend = F, legend_pos = "topleft", main_header = "", 
                     pos_main = 0.0){

  
  col_mea <- "steelblue4"
  col_obs <- "darkred"
  col_clm <- "grey25"
  
  doy_df <- data.frame(doy_mea = range(doy_mea, na.rm = T),
                       doy_obs = range(doy_obs, na.rm = T),
                       doy_cm1 = range(doy_cm1, na.rm = T),
                       doy_cm2 = range(doy_cm2, na.rm = T),
                       doy_cm3 = range(doy_cm3, na.rm = T),
                       doy_cm4 = range(doy_cm4, na.rm = T),
                       doy_cm5 = range(doy_cm5, na.rm = T))
  
  boxplot(doy_df, boxfill = NA, border = NA, axes = F, ylim = ylims)
  axis(2, mgp=c(3, 0.55, 0), tck = -0.017, cex.axis = 2.5)
  mtext(y_lab, side = 2, line = 3.1, cex = 1.8)
  grid(nx = 0, ny = 6, col = "grey55", lwd = 0.5)
  if(do_legend){
    legend(legend_pos, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
           col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.3,
           box.lwd = 0.0, box.col = "black", bg = "white")
  }
  box()
  
  boxplot(doy_mea, ylim = ylims, col = col_mea, axes = F, xaxt = "n", add = TRUE, 
          at = 1, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_obs, ylim = ylims, col = col_obs, axes = F, xaxt = "n", add = TRUE, 
          at = 2, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_cm1, ylim = ylims, col = col_clm, axes = F, xaxt = "n", add = TRUE, 
          at = 3, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_cm2, ylim = ylims, col = col_clm, axes = F, xaxt = "n", add = TRUE, 
          at = 4, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_cm3, ylim = ylims, col = col_clm, axes = F, xaxt = "n", add = TRUE, 
          at = 5, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_cm4, ylim = ylims, col = col_clm, axes = F, xaxt = "n", add = TRUE, 
          at = 6, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(doy_cm5, ylim = ylims, col = col_clm, axes = F, xaxt = "n", add = TRUE, 
          at = 7, boxwex = 1.5, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  
  
  mtext(main_header, side = 3, line = 0.5, adj = pos_main, cex = 2.2)
  
}

plot_doy(doy_mea = disc_mea_doy[, 1], doy_obs = disc_obs_doy[, 1], doy_cm1 = disc_cm1_doy[, 1], 
         doy_cm2 = disc_cm2_doy[, 1], doy_cm3 = disc_cm3_doy[, 1], 
         doy_cm4 = disc_cm4_doy[, 1], doy_cm5 = disc_cm5_doy[, 1],
         main_header = "a) Lobith")

plot_doy(doy_mea = disc_mea_doy[, 2], doy_obs = disc_obs_doy[, 2], doy_cm1 = disc_cm1_doy[, 2], 
         doy_cm2 = disc_cm2_doy[, 2], doy_cm3 = disc_cm3_doy[, 2], 
         doy_cm4 = disc_cm4_doy[, 2], doy_cm5 = disc_cm5_doy[, 2],
         main_header = "b) Cologne")

plot_doy(doy_mea = disc_mea_doy[, 3], doy_obs = disc_obs_doy[, 3], doy_cm1 = disc_cm1_doy[, 3], 
         doy_cm2 = disc_cm2_doy[, 3], doy_cm3 = disc_cm3_doy[, 3], 
         doy_cm4 = disc_cm4_doy[, 3], doy_cm5 = disc_cm5_doy[, 3],
         main_header = "c) Cochem")

plot_doy(doy_mea = disc_mea_doy[, 4], doy_obs = disc_obs_doy[, 4], doy_cm1 = disc_cm1_doy[, 4], 
         doy_cm2 = disc_cm2_doy[, 4], doy_cm3 = disc_cm3_doy[, 4], 
         doy_cm4 = disc_cm4_doy[, 4], doy_cm5 = disc_cm5_doy[, 4],
         main_header = "d) Kaub")

plot_doy(doy_mea = disc_mea_doy[, 5], doy_obs = disc_obs_doy[, 5], doy_cm1 = disc_cm1_doy[, 5], 
         doy_cm2 = disc_cm2_doy[, 5], doy_cm3 = disc_cm3_doy[, 5], 
         doy_cm4 = disc_cm4_doy[, 5], doy_cm5 = disc_cm5_doy[, 5],
         main_header = "e) Wuerzburg")

plot_doy(doy_mea = disc_mea_doy[, 6], doy_obs = disc_obs_doy[, 6], doy_cm1 = disc_cm1_doy[, 6], 
         doy_cm2 = disc_cm2_doy[, 6], doy_cm3 = disc_cm3_doy[, 6], 
         doy_cm4 = disc_cm4_doy[, 6], doy_cm5 = disc_cm5_doy[, 6],
         main_header = "f) Worms")

plot_doy(doy_mea = disc_mea_doy[, 7], doy_obs = disc_obs_doy[, 7], doy_cm1 = disc_cm1_doy[, 7], 
         doy_cm2 = disc_cm2_doy[, 7], doy_cm3 = disc_cm3_doy[, 7], 
         doy_cm4 = disc_cm4_doy[, 7], doy_cm5 = disc_cm5_doy[, 7],
         main_header = "g) Rockenau")

plot_doy(doy_mea = disc_mea_doy[, 8], doy_obs = disc_obs_doy[, 8], doy_cm1 = disc_cm1_doy[, 8], 
         doy_cm2 = disc_cm2_doy[, 8], doy_cm3 = disc_cm3_doy[, 8], 
         doy_cm4 = disc_cm4_doy[, 8], doy_cm5 = disc_cm5_doy[, 8],
         main_header = "h) Speyer")

plot_doy(doy_mea = disc_mea_doy[, 9], doy_obs = disc_obs_doy[, 9], doy_cm1 = disc_cm1_doy[, 9], 
         doy_cm2 = disc_cm2_doy[, 9], doy_cm3 = disc_cm3_doy[, 9], 
         doy_cm4 = disc_cm4_doy[, 9], doy_cm5 = disc_cm5_doy[, 9],
         main_header = "i) Basel")

#Label y axis
par(mar = c(0, 0, 0, 0))
plot(1:10, 1:10, type = "n", ylab = "", xlab = "", axes = F)
mtext("Timing annual maxima [DOY]", side = 2, line = -3.5, cex = 2.5)

dev.off()



#monthly_quantiles----

perc_month <- function(data_in){
  
  disc_day <- ord_day(data_in = data_in,
                      date = date_simu,
                      start_y = 1951,
                      end_y = 2000,
                      break_day = 274)
  
  jan_cols <- 1:31
  feb_cols <- 32:59
  mar_cols <- 60:90
  apr_cols <- 91:120
  may_cols <- 121:151
  jun_cols <- 152:181
  jul_cols <- 182:212
  aug_cols <- 213:243
  sep_cols <- 244:273
  oct_cols <- 274:304
  nov_cols <- 305:334
  dec_cols <- 335:365
  
  month_cols <- list(jan_cols, feb_cols, mar_cols, apr_cols, may_cols, jun_cols, jul_cols, aug_cols, sep_cols, oct_cols, nov_cols, dec_cols)
  perc_sel <- 0.90
  
  qmon <- rep(NA, 12)
  
  for(i in 1:12){
    qmon[i] <- stats::quantile(disc_day[ , month_cols[[i]]], probs = perc_sel, 
                                        type = 8, na.rm = T)
  }
  
  return(qmon)
  
}

qmon_cm1 <- apply(disc_cm1, 2, perc_month)
qmon_cm2 <- apply(disc_cm2, 2, perc_month)
qmon_cm3 <- apply(disc_cm3, 2, perc_month)
qmon_cm4 <- apply(disc_cm4, 2, perc_month)
qmon_cm5 <- apply(disc_cm5, 2, perc_month)
qmon_eob <- apply(disc_obs, 2, perc_month)
qmon_obs <- cbind(perc_month(disc_lobi$value), perc_month(disc_koel$value), perc_month(disc_coch$value),
                  perc_month(disc_kaub$value), perc_month(disc_wuer$value), perc_month(disc_worm$value),
                  perc_month(disc_rock$value), perc_month(disc_spey$value), perc_month(disc_base$value))

pdf(paste0(bas_dir,"res_figs/eval_hist_qua_raw.pdf"), width = 16, height = 8)

layout(matrix(c(10, 1, 2, 3,
                10, 4, 5, 6,
                10, 7, 8, 9),
              3, 4, byrow = T), widths=c(0.08, 1, 1, 1), heights=c())

par(family = "serif")
par(mar = c(2.5, 3.5, 3.0, 1.0))

plot_month_quan <- function(col_sel, main = ""){
  
  ylims <- range(qmon_obs[, col_sel], qmon_eob[, col_sel], qmon_cm1[, col_sel], qmon_cm2[, col_sel], 
                 qmon_cm3[, col_sel], qmon_cm4[, col_sel], qmon_cm5[, col_sel])
  
  x_labs <- c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S")
  x_ats <- seq(1, 12, 1)
  x_tics <- seq(0.5, 12.5, 1)
  line_lwd <- 2.0
  point_cex <- 1.4
  
  col_gcm <- alpha("black", alpha = 0.5)
  col_meas <- alpha("steelblue4", alpha = 0.5)
  col_eobs <- alpha("darkred", alpha = 0.5)
  
  plot(qmon_obs[, col_sel], type = "n", ylim = ylims, axes = F, ylab = "", 
       xlab = "", xlim = c(0.5, 12.5))
  grid(nx = 0, ny = 4, col = "grey75", lty = "dashed", lwd = 0.7)
  lines(qmon_cm1[, col_sel], col = col_gcm, lwd = line_lwd)
  lines(qmon_cm2[, col_sel], col = col_gcm, lwd = line_lwd)
  lines(qmon_cm3[, col_sel], col = col_gcm, lwd = line_lwd)
  lines(qmon_cm4[, col_sel], col = col_gcm, lwd = line_lwd)
  lines(qmon_cm5[, col_sel], col = col_gcm, lwd = line_lwd)
  lines(qmon_obs[, col_sel], col = col_meas, lwd = line_lwd)
  lines(qmon_eob[, col_sel], col = col_eobs, lwd = line_lwd)
  points(qmon_cm1[, col_sel], bg = "grey25", col = "grey25", pch = 21, cex = point_cex)
  points(qmon_cm2[, col_sel], bg = "grey25", col = "grey25", pch = 22, cex = point_cex)
  points(qmon_cm3[, col_sel], bg = "grey25", col = "grey25", pch = 23, cex = point_cex)
  points(qmon_cm4[, col_sel], bg = "grey25", col = "grey25", pch = 24, cex = point_cex)
  points(qmon_cm5[, col_sel], bg = "grey25", col = "grey25", pch = 25, cex = point_cex)
  points(qmon_obs[, col_sel], col = "steelblue4", pch = 19, cex = point_cex)
  points(qmon_eob[, col_sel], col = "darkred", pch = 19, cex = point_cex)
  axis(1, labels = x_labs, at = x_ats, mgp=c(3, 0.69, 0), tck = -0.0, cex.axis = 2.0)
  axis(1, labels = F, at = x_tics, tck = -0.050)
  axis(2, mgp=c(3, 0.21, 0), tck = -0.015, cex.axis = 2.5)
  mtext(main, side = 3, adj = 0, cex = 2.2, line = 0.5)
  box()
  
}

plot_month_quan(1, "a) Lobith")

plot_month_quan(2, "b) Cologne")

plot_month_quan(3, "c) Cochem")

plot_month_quan(4, "d) Kaub")

plot_month_quan(5, "e) Wuerzburg")

plot_month_quan(6, "f) Worms")

plot_month_quan(7, "g) Rockenau")

plot_month_quan(8, "h) Speyer")

plot_month_quan(9, "i) Basel")

#Label y axis
par(mar = c(0, 0, 0, 0))
plot(1:10, 1:10, type = "n", ylab = "", xlab = "", axes = F)
mtext(expression(paste("90 % Quantile discharge", " [m"^"3", "s"^"-1","]")), 
      side = 2, line = -3.6, cex = 2.4)

dev.off()

