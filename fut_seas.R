###

#mHM simulations Part II: mHM simulations future scenarios RCP2.6, RCP6.0 and RCP8.5
#Erwin Rottler, UniversityI of Potsdam, Summer 2020

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales, ncdf4,
               ncdf4.helpers, sp, raster, viridis, meltimr, POT, readr, hydroGOF,
               CoinCalc, seas, tictoc, lubridate, plyr)

#set directories
bas_dir <- "U:/rhine_fut/R/"
run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

#load functions
source(paste0(bas_dir, "mhm_rhine/functs.R"))

# stopCluster(my_clust)
# 
# n_cores <- 5 #number of cores used for parallel computing
# 
# #Make cluster for parallel computing
# my_clust <- makeCluster(n_cores)
# clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
# clusterExport(my_clust, "index_col_base")
# registerDoParallel(my_clust)

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#read table warming levels
warming_periods <- read.table("D:/nrc_user/rottler/mhm_run/6435060/input/meteo/periods_warming.txt",
                              header = T, stringsAsFactors = F)

#paths to mHM nc-files with fluxes and states
nc_file_paths <- 
  c(paste0(run_dir, "output/GFDL-ESM2M/historical/"),
    paste0(run_dir, "output/HadGEM2-ES/historical/"),
    paste0(run_dir, "output/IPSL-CM5A-LR/historical/"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/historical/"),
    paste0(run_dir, "output/NorESM1-M/historical/"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp2p6/"),
    paste0(run_dir, "output/HadGEM2-ES/rcp2p6/"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp2p6/"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp2p6/"),
    paste0(run_dir, "output/NorESM1-M/rcp2p6/"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp6p0/"),
    paste0(run_dir, "output/HadGEM2-ES/rcp6p0/"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp6p0/"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp6p0/"),
    paste0(run_dir, "output/NorESM1-M/rcp6p0/"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp8p5/"),
    paste0(run_dir, "output/HadGEM2-ES/rcp8p5/"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp8p5/"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp8p5/"),
    paste0(run_dir, "output/NorESM1-M/rcp8p5/"))

#sel_gauges----

#get coordinates for selected gauges
koel_file <- paste0(grdc_dir, "6335060_Q_Day.Cmd.txt")
coch_file <- paste0(grdc_dir, "6336050_Q_Day.Cmd.txt")
base_file <- paste0(grdc_dir, "6935051_Q_Day.Cmd.txt")

file_paths <- c(koel_file, coch_file, base_file)

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

#Get grid cells of selected statins
nc_disc_file_obs <- paste0(run_dir, "output/EOBS/output/mRM_Fluxes_States.nc")
nc_disc_obs <- ncdf4::nc_open(nc_disc_file_obs)

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

rows_sel_gaugs[2] <- rows_sel_gaugs[2]+1 #gauge cochem one row lower


#get_disc----

#paths to mRM nc-files with routed discharge
nc_file_disc_paths <- 
  c(paste0(run_dir, "output/GFDL-ESM2M/historical/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/historical/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/historical/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/historical/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/historical/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp2p6/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp2p6/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp2p6/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp2p6/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp2p6/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp6p0/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp6p0/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp6p0/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp6p0/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp6p0/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp8p5/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp8p5/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp8p5/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp8p5/output/mRM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp8p5/output/mRM_Fluxes_States.nc"))

#Function to extract time series from nc-files
disc_from_nc <- function(gcm_model, delta_t, rcp){
  
  #select nc_file
  nc_path_sel <- nc_file_disc_paths[which(grepl(gcm_model, nc_file_disc_paths) &
                                            grepl(rcp, nc_file_disc_paths))]
  nc_file_sel <- nc_open(nc_path_sel)
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  count_date <- length(date)
  
  simu_koel <- ncvar_get(nc_file_sel, start = c(rows_sel_gaugs[1], cols_sel_gaugs[1], 1), 
                         count = c(1, 1, count_date), varid = "Qrouted")
  simu_coch <- ncvar_get(nc_file_sel, start = c(rows_sel_gaugs[2], cols_sel_gaugs[2], 1), 
                         count = c(1, 1, count_date), varid = "Qrouted")
  simu_base <- ncvar_get(nc_file_sel, start = c(rows_sel_gaugs[3], cols_sel_gaugs[3], 1), 
                         count = c(1, 1, count_date), varid = "Qrouted")
  
  disc_sim_full <- cbind(simu_koel, simu_coch, simu_base)
  
  disc_sim <- disc_sim_full[which(date %in% date_sel), ]
  
  #add date column
  date_out <- date[which(date %in% date_sel)]
  
  disc_out <- data.frame(date = date_out,
                         koel = disc_sim[, 1],
                         coch = disc_sim[, 2],
                         base = disc_sim[, 3])
  
  return(disc_out)
  
}

#historical
tic()
disc_hist_1 <- disc_from_nc("GFDL-ESM2M",     "historical", "historical")
disc_hist_2 <- disc_from_nc("HadGEM2-ES",     "historical", "historical")
disc_hist_3 <- disc_from_nc("IPSL-CM5A-LR",   "historical", "historical")
disc_hist_4 <- disc_from_nc("MIROC-ESM-CHEM", "historical", "historical")
disc_hist_5 <- disc_from_nc("NorESM1-M",      "historical", "historical")
toc()

#1.5K warming level
tic()
disc_1p5K_1  <- disc_from_nc("HadGEM2-ES",     "1p5", "2p6")
disc_1p5K_2  <- disc_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
disc_1p5K_3  <- disc_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
disc_1p5K_4  <- disc_from_nc("NorESM1-M",      "1p5", "2p6")
disc_1p5K_5  <- disc_from_nc("GFDL-ESM2M",     "1p5", "6p0")
disc_1p5K_6  <- disc_from_nc("HadGEM2-ES",     "1p5", "6p0")
disc_1p5K_7  <- disc_from_nc("IPSL-CM5A-LR",   "1p5", "6p0")
disc_1p5K_8  <- disc_from_nc("MIROC-ESM-CHEM", "1p5", "6p0")
disc_1p5K_9  <- disc_from_nc("NorESM1-M",      "1p5", "6p0")
disc_1p5K_10 <- disc_from_nc("GFDL-ESM2M",     "1p5", "8p5")
disc_1p5K_11 <- disc_from_nc("HadGEM2-ES",     "1p5", "8p5")
disc_1p5K_12 <- disc_from_nc("IPSL-CM5A-LR",   "1p5", "8p5")
disc_1p5K_13 <- disc_from_nc("MIROC-ESM-CHEM", "1p5", "8p5")
disc_1p5K_14 <- disc_from_nc("NorESM1-M",      "1p5", "8p5")
toc()

#2K warming level
tic()
disc_2p0K_1  <- disc_from_nc("HadGEM2-ES",     "2p0", "2p6")
disc_2p0K_2  <- disc_from_nc("IPSL-CM5A-LR",   "2p0", "2p6")
disc_2p0K_3  <- disc_from_nc("MIROC-ESM-CHEM", "2p0", "2p6")
disc_2p0K_4  <- disc_from_nc("GFDL-ESM2M",     "2p0", "6p0")
disc_2p0K_5  <- disc_from_nc("HadGEM2-ES",     "2p0", "6p0")
disc_2p0K_6  <- disc_from_nc("IPSL-CM5A-LR",   "2p0", "6p0")
disc_2p0K_7  <- disc_from_nc("MIROC-ESM-CHEM", "2p0", "6p0")
disc_2p0K_8  <- disc_from_nc("NorESM1-M",      "2p0", "6p0")
disc_2p0K_9  <- disc_from_nc("GFDL-ESM2M",     "2p0", "8p5")
disc_2p0K_10 <- disc_from_nc("HadGEM2-ES",     "2p0", "8p5")
disc_2p0K_11 <- disc_from_nc("IPSL-CM5A-LR",   "2p0", "8p5")
disc_2p0K_12 <- disc_from_nc("MIROC-ESM-CHEM", "2p0", "8p5")
disc_2p0K_13 <- disc_from_nc("NorESM1-M",      "2p0", "8p5")
toc()

#3K warming level
tic()
disc_3p0K_1 <- disc_from_nc("HadGEM2-ES",     "3p0", "6p0")
disc_3p0K_2 <- disc_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
disc_3p0K_3 <- disc_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
disc_3p0K_4 <- disc_from_nc("GFDL-ESM2M",     "3p0", "8p5")
disc_3p0K_5 <- disc_from_nc("HadGEM2-ES",     "3p0", "8p5")
disc_3p0K_6 <- disc_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
disc_3p0K_7 <- disc_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
disc_3p0K_8 <- disc_from_nc("NorESM1-M",      "3p0", "8p5")
toc()

#get_fluxes----

# dem <- raster(paste0(run_dir, "input/morph/dem.asc"), crs = epsg3035)
dem = raster("U:/rhine_snow/data/eu_dem/processed/eu_dem_1000.tif")

#Get basin .shp
basin_coch_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/cochem_catch.shp")
basin_base_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/basel_catch.shp")
basin_koel_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/koeln_catch.shp")
basin_coch <- spTransform(basin_coch_raw, CRS = crswgs84)
basin_base <- spTransform(basin_base_raw, CRS = crswgs84)
basin_koel <- spTransform(basin_koel_raw, CRS = crswgs84)

#get cells in basin
nc_dummy <- nc_open(paste0(run_dir, "output/GFDL-ESM2M/historical/output/mRM_Fluxes_States.nc"))
lon <- ncdf4::ncvar_get(nc_dummy, varid = "lon")
lat <- ncdf4::ncvar_get(nc_dummy, varid = "lat")

#spatial grid points from lat/lon
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)

#grid points inside watersheds
inside_base <- !is.na(sp::over(grid_points_cube_84, as(basin_base, "SpatialPolygons")))
inside_coch <- !is.na(sp::over(grid_points_cube_84, as(basin_coch, "SpatialPolygons")))
inside_koel <- !is.na(sp::over(grid_points_cube_84, as(basin_koel, "SpatialPolygons")))

grid_points_base <- grid_points_cube_84[which(inside_base == T)]
grid_points_coch <- grid_points_cube_84[which(inside_coch == T)]
grid_points_koel <- grid_points_cube_84[which(inside_koel == T)]

#Select cells for Basel/Cochem watershed
lat_in_base <- grid_points_base@coords[, 2]
lat_in_coch <- grid_points_coch@coords[, 2]
lat_in_koel <- grid_points_koel@coords[, 2]

my_get_cube_col <- function(val_in, coor_in = lat){
  
  get_cube_index_col(val_in = val_in, coor_in = coor_in)
  
}
my_get_cube_row <- function(val_in, coor_in = lat){
  
  get_cube_index_row(val_in = val_in, coor_in = coor_in)
  
}

#get index in cube from points inside sub-basins
index_col_base <- sapply(lat_in_base, my_get_cube_col)
index_row_base <- sapply(lat_in_base, my_get_cube_row)
index_col_coch <- sapply(lat_in_coch, my_get_cube_col)
index_row_coch <- sapply(lat_in_coch, my_get_cube_row)
index_col_koel <- sapply(lat_in_koel, my_get_cube_col)
index_row_koel <- sapply(lat_in_koel, my_get_cube_row)

#elevations
grid_points_base_dem <- spTransform(grid_points_base, crs(dem))
grid_points_coch_dem <- spTransform(grid_points_coch, crs(dem))
grid_points_koel_dem <- spTransform(grid_points_koel, crs(dem))

n_cores <- 5 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

elevs_base <- foreach(i = 1:length(grid_points_base_dem), .combine = 'c') %dopar% {
  
  elev_buff(grid_points_base_dem[i], 2500, dem_in = dem)
  
}

elevs_coch <- foreach(i = 1:length(grid_points_coch_dem), .combine = 'c') %dopar% {
  
  elev_buff(grid_points_coch_dem[i], 2500, dem_in = dem)
  
}

elevs_koel <- foreach(i = 1:length(grid_points_koel_dem), .combine = 'c') %dopar% {
  
  elev_buff(grid_points_koel_dem[i], 2500, dem_in = dem)
  
}

stopCluster(my_clust)

#Function to extract time series from nc-files
flux_from_nc <- function(gcm_model, delta_t, rcp, ncores = 5){
  
  #select nc_file
  nc_path_sel <- nc_file_paths[which(grepl(gcm_model, nc_file_paths) &
                                          grepl(rcp, nc_file_paths))]
  nc_file_sel <- nc_open(paste0(nc_path_sel, "output/mHM_Fluxes_States.nc"))
   
  #path to precipitation input
  # nc_path_sel_prec <- gsub("output", "input/meteo", nc_path_sel)
  # pr_file_sel <- nc_open(paste0(nc_path_sel, "input/pre.nc"))
                                                
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  
  #if simulation time frame does not entirely cover warming period
  if(date_sel[1] > date[1]){
    sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
    count_date <- length(date_sel)
  }else{
    sta_date_ind <- 1
    count_date <- length(date_sel) - which(format(date_sel) == date[1]) + 1
  }
  
  #snowpack
  snow_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
  #effective precipitation
  pef_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lon), count_date), varid = "preEffect")
  
  n_cores <- ncores #number of cores used for parallel computing
  
  #Make cluster for parallel computing
  my_clust <- makeCluster(n_cores)
  clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
  clusterExport(my_clust, c("index_col_base", "index_row_base", "index_col_coch", "index_row_coch",
                            "index_col_koel", "index_row_koel", "elevs_base", "elevs_coch", "elevs_koel"))
  registerDoParallel(my_clust)
  
  sno_base <- foreach(i = 1:length(index_col_base), .combine = 'cbind') %dopar% {
    snow_cube[index_col_base[i], index_row_base[i], ]
  }
  sno_coch <- foreach(i = 1:length(index_col_coch), .combine = 'cbind') %dopar% {
    snow_cube[index_col_coch[i], index_row_coch[i], ]
  }
  sno_koel <- foreach(i = 1:length(index_col_koel), .combine = 'cbind') %dopar% {
    snow_cube[index_col_koel[i], index_row_koel[i], ]
  }
  
  epn_base <- foreach(i = 1:length(index_col_base), .combine = 'cbind') %dopar% {
    pef_cube[index_col_base[i], index_row_base[i], ]
  }
  epn_coch <- foreach(i = 1:length(index_col_coch), .combine = 'cbind') %dopar% {
    pef_cube[index_col_coch[i], index_row_coch[i], ]
  }
  epn_koel <- foreach(i = 1:length(index_col_koel), .combine = 'cbind') %dopar% {
    pef_cube[index_col_koel[i], index_row_koel[i], ]
  }
  
  stopCluster(my_clust)
  
  #Values on basin scale
  base_sd_mea <- apply(sno_base, 1, mea_na)
  base_sd_mea_dif <- c(NA, diff(base_sd_mea))
  base_sd_mea_dif[which(base_sd_mea_dif > 0)] <- NA
  base_me_mea <- base_sd_mea_dif * -1 #melt positive values
  base_ep_mea <- apply(epn_base, 1, mea_na)
  
  coch_sd_mea <- apply(sno_coch, 1, mea_na)
  coch_sd_mea_dif <- c(NA, diff(coch_sd_mea))
  coch_sd_mea_dif[which(coch_sd_mea_dif > 0)] <- NA
  coch_me_mea <- coch_sd_mea_dif * -1 #melt positive values
  coch_ep_mea <- apply(epn_coch, 1, mea_na)
  
  koel_sd_mea <- apply(sno_koel, 1, mea_na)
  koel_sd_mea_dif <- c(NA, diff(koel_sd_mea))
  koel_sd_mea_dif[which(koel_sd_mea_dif > 0)] <- NA
  koel_me_mea <- koel_sd_mea_dif * -1 #melt positive values
  koel_ep_mea <- apply(epn_koel, 1, mea_na)
  
  base_melt_ma <- rollapply(data = base_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
  coch_melt_ma <- rollapply(data = coch_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
  koel_melt_ma <- rollapply(data = koel_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
  
  #Moving window liquid precipitation
  base_lpre <- base_ep_mea - base_me_mea
  coch_lpre <- coch_ep_mea - coch_me_mea
  koel_lpre <- koel_ep_mea - koel_me_mea
  base_lpre_ma <- rollapply(data = base_lpre, width = 5,
                            FUN = sum_na, align = "right", fill = NA)
  coch_lpre_ma <- rollapply(data = coch_lpre, width = 5,
                            FUN = sum_na, align = "right", fill = NA)
  koel_lpre_ma <- rollapply(data = koel_lpre, width = 5,
                            FUN = sum_na, align = "right", fill = NA)
  
  #fraction snowmelt contribution
  base_frac <- base_melt_ma / (base_melt_ma + base_lpre_ma)
  coch_frac <- coch_melt_ma / (coch_melt_ma + coch_lpre_ma)
  koel_frac <- koel_melt_ma / (koel_melt_ma + koel_lpre_ma)
  
  date_out <- date[which(date %in% date_sel)]
  
  #get mean melt elevation
  f_melt_calc <- function(snow_depth){
    
    snow_depth_dif <- c(NA, diff(snow_depth))
    snow_depth_dif[which(snow_depth_dif > 0)] <- NA
    melt <- snow_depth_dif * -1 #melt positive values
    
    return(melt)
    
  }
  
  melt_base <- apply(sno_base, 2, f_melt_calc)
  melt_coch <- apply(sno_coch, 2, f_melt_calc)
  melt_koel <- apply(sno_koel, 2, f_melt_calc)
  
  f_mea_ele_base <- function(melt_in){
    sum_na((melt_in / sum_na(melt_in)) * elevs_base)
  }
  f_mea_ele_coch <- function(melt_in){
    sum_na((melt_in / sum_na(melt_in)) * elevs_coch)
  }
  f_mea_ele_koel <- function(melt_in){
    sum_na((melt_in / sum_na(melt_in)) * elevs_koel)
  }
  
  mel_ele_base <- apply(melt_base, 1, f_mea_ele_base)
  mel_ele_coch <- apply(melt_coch, 1, f_mea_ele_coch)
  mel_ele_koel <- apply(melt_koel, 1, f_mea_ele_koel)
  
  flux_out <- data.frame(date = date_out,
                         base_melt_ma = base_melt_ma,
                         coch_melt_ma = coch_melt_ma,
                         koel_melt_ma = koel_melt_ma,
                         base_lpre_ma = base_lpre_ma,
                         coch_lpre_ma = coch_lpre_ma,
                         koel_lpre_ma = koel_lpre_ma,
                         base_frac = base_frac,
                         coch_frac = coch_frac,
                         koel_frac = koel_frac,
                         base_mel_ele = mel_ele_base,
                         coch_mel_ele = mel_ele_coch,
                         koel_mel_ele = mel_ele_koel)  
  
  return(flux_out)
  
}
prec_from_nc <- function(gcm_model, delta_t, rcp, ncores = 5){
  
  #select nc_file
  nc_path_sel <- nc_file_paths[which(grepl(gcm_model, nc_file_paths) &
                                       grepl(rcp, nc_file_paths))]
  
  nc_path_sel_meteo <- gsub("output", "input/meteo", nc_path_sel)
  
  nc_file_pre <- nc_open(paste0(nc_path_sel_meteo, "pre.nc"))
  nc_file_tem <- nc_open(paste0(nc_path_sel_meteo, "tavg.nc"))
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_pre, time.dim.name = "time")))
  
  #if simulation time frame does not entirely cover warming period
  if(date_sel[1] > date[1]){
    sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
    count_date <- length(date_sel)
  }else{
    sta_date_ind <- 1
    count_date <- length(date_sel) - which(format(date_sel) == date[1]) + 1
  }
  
  #precipitation
  prec_cube <- ncvar_get(nc_file_pre, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "pre")
  #temperature
  temp_cube <- ncvar_get(nc_file_tem, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lon), count_date), varid = "tavg")
  
  n_cores <- ncores #number of cores used for parallel computing
  
  #Make cluster for parallel computing
  my_clust <- makeCluster(n_cores)
  clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
  clusterExport(my_clust, c("index_col_base", "index_row_base", "index_col_coch", "index_row_coch",
                            "index_col_koel", "index_row_koel", "elevs_base", "elevs_coch", "elevs_koel"))
  registerDoParallel(my_clust)
  
  pre_base <- foreach(i = 1:length(index_col_base), .combine = 'cbind') %dopar% {
    prec_cube[index_col_base[i], index_row_base[i], ]
  }
  pre_coch <- foreach(i = 1:length(index_col_coch), .combine = 'cbind') %dopar% {
    prec_cube[index_col_coch[i], index_row_coch[i], ]
  }
  pre_koel <- foreach(i = 1:length(index_col_koel), .combine = 'cbind') %dopar% {
    prec_cube[index_col_koel[i], index_row_koel[i], ]
  }
  
  tem_base <- foreach(i = 1:length(index_col_base), .combine = 'cbind') %dopar% {
    temp_cube[index_col_base[i], index_row_base[i], ]
  }
  tem_coch <- foreach(i = 1:length(index_col_coch), .combine = 'cbind') %dopar% {
    temp_cube[index_col_coch[i], index_row_coch[i], ]
  }
  tem_koel <- foreach(i = 1:length(index_col_koel), .combine = 'cbind') %dopar% {
    temp_cube[index_col_koel[i], index_row_koel[i], ]
  }
  
  stopCluster(my_clust)
  
  #total precipitation
  base_pre_mea <- apply(pre_base, 1, mea_na)
  coch_pre_mea <- apply(pre_coch, 1, mea_na)
  koel_pre_mea <- apply(pre_koel, 1, mea_na)
  
  #liquid precipitation
  temp_thresh <- 0.965483531537 #determined at calibration
  
  for(i in 1:ncol(pre_base)){
    pre_base[which(tem_base[, i] < temp_thresh), i] <- 0
  }
  
  for(i in 1:ncol(pre_coch)){
    pre_coch[which(tem_coch[, i] < temp_thresh), i] <- 0
  }
  
  for(i in 1:ncol(pre_koel)){
    pre_koel[which(tem_koel[, i] < temp_thresh), i] <- 0
  }
  
  base_lpr_mea <- apply(pre_base, 1, mea_na)
  coch_lpr_mea <- apply(pre_coch, 1, mea_na)
  koel_lpr_mea <- apply(pre_koel, 1, mea_na)
  
  #protective effect
  base_pro_eff <- 1 - (base_lpr_mea / base_pre_mea)
  coch_pro_eff <- 1 - (coch_lpr_mea / coch_pre_mea)
  koel_pro_eff <- 1 - (koel_lpr_mea / koel_pre_mea)
  
  date_out <- date[which(date %in% date_sel)]
  
  prec_out <- data.frame(date = date_out,
                         base_pre_tot = base_pre_mea,
                         coch_pre_tot = coch_pre_mea,
                         koel_pre_tot = koel_pre_mea,
                         base_pre_liq = base_lpr_mea,
                         coch_pre_liq = coch_lpr_mea,
                         koel_pre_liq = koel_lpr_mea,
                         base_pro_eff = base_pro_eff,
                         coch_pro_eff = coch_pro_eff,
                         koel_pro_eff = koel_pro_eff)
  
}

#historical
tic()
flux_hist_1 <- flux_from_nc("GFDL-ESM2M",     "historical", "historical")
flux_hist_2 <- flux_from_nc("HadGEM2-ES",     "historical", "historical")
flux_hist_3 <- flux_from_nc("IPSL-CM5A-LR",   "historical", "historical")
flux_hist_4 <- flux_from_nc("MIROC-ESM-CHEM", "historical", "historical")
flux_hist_5 <- flux_from_nc("NorESM1-M",      "historical", "historical")
toc()

tic()
prec_hist_1 <- prec_from_nc("GFDL-ESM2M",     "historical", "historical")
prec_hist_2 <- prec_from_nc("HadGEM2-ES",     "historical", "historical")
prec_hist_3 <- prec_from_nc("IPSL-CM5A-LR",   "historical", "historical")
prec_hist_4 <- prec_from_nc("MIROC-ESM-CHEM", "historical", "historical")
prec_hist_5 <- prec_from_nc("NorESM1-M",      "historical", "historical")
toc()

#1.5K warming level
tic()
flux_1p5K_1  <- flux_from_nc("HadGEM2-ES",     "1p5", "2p6")
flux_1p5K_2  <- flux_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
flux_1p5K_3  <- flux_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
flux_1p5K_4  <- flux_from_nc("NorESM1-M",      "1p5", "2p6")
flux_1p5K_5  <- flux_from_nc("GFDL-ESM2M",     "1p5", "6p0")
flux_1p5K_6  <- flux_from_nc("HadGEM2-ES",     "1p5", "6p0")
flux_1p5K_7  <- flux_from_nc("IPSL-CM5A-LR",   "1p5", "6p0")
flux_1p5K_8  <- flux_from_nc("MIROC-ESM-CHEM", "1p5", "6p0")
flux_1p5K_9  <- flux_from_nc("NorESM1-M",      "1p5", "6p0")
flux_1p5K_10 <- flux_from_nc("GFDL-ESM2M",     "1p5", "8p5")
flux_1p5K_11 <- flux_from_nc("HadGEM2-ES",     "1p5", "8p5")
flux_1p5K_12 <- flux_from_nc("IPSL-CM5A-LR",   "1p5", "8p5")
flux_1p5K_13 <- flux_from_nc("MIROC-ESM-CHEM", "1p5", "8p5")
flux_1p5K_14 <- flux_from_nc("NorESM1-M",      "1p5", "8p5")
toc()

tic()
prec_1p5K_1  <- prec_from_nc("HadGEM2-ES",     "1p5", "2p6")
prec_1p5K_2  <- prec_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
prec_1p5K_3  <- prec_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
prec_1p5K_4  <- prec_from_nc("NorESM1-M",      "1p5", "2p6")
prec_1p5K_5  <- prec_from_nc("GFDL-ESM2M",     "1p5", "6p0")
prec_1p5K_6  <- prec_from_nc("HadGEM2-ES",     "1p5", "6p0")
prec_1p5K_7  <- prec_from_nc("IPSL-CM5A-LR",   "1p5", "6p0")
prec_1p5K_8  <- prec_from_nc("MIROC-ESM-CHEM", "1p5", "6p0")
prec_1p5K_9  <- prec_from_nc("NorESM1-M",      "1p5", "6p0")
prec_1p5K_10 <- prec_from_nc("GFDL-ESM2M",     "1p5", "8p5")
prec_1p5K_11 <- prec_from_nc("HadGEM2-ES",     "1p5", "8p5")
prec_1p5K_12 <- prec_from_nc("IPSL-CM5A-LR",   "1p5", "8p5")
prec_1p5K_13 <- prec_from_nc("MIROC-ESM-CHEM", "1p5", "8p5")
prec_1p5K_14 <- prec_from_nc("NorESM1-M",      "1p5", "8p5")
toc()

#2K warming level
tic()
flux_2p0K_1  <- flux_from_nc("HadGEM2-ES",     "2p0", "2p6")
flux_2p0K_2  <- flux_from_nc("IPSL-CM5A-LR",   "2p0", "2p6")
flux_2p0K_3  <- flux_from_nc("MIROC-ESM-CHEM", "2p0", "2p6")
flux_2p0K_4  <- flux_from_nc("GFDL-ESM2M",     "2p0", "6p0")
flux_2p0K_5  <- flux_from_nc("HadGEM2-ES",     "2p0", "6p0")
flux_2p0K_6  <- flux_from_nc("IPSL-CM5A-LR",   "2p0", "6p0")
flux_2p0K_7  <- flux_from_nc("MIROC-ESM-CHEM", "2p0", "6p0")
flux_2p0K_8  <- flux_from_nc("NorESM1-M",      "2p0", "6p0")
flux_2p0K_9  <- flux_from_nc("GFDL-ESM2M",     "2p0", "8p5")
flux_2p0K_10 <- flux_from_nc("HadGEM2-ES",     "2p0", "8p5")
flux_2p0K_11 <- flux_from_nc("IPSL-CM5A-LR",   "2p0", "8p5")
flux_2p0K_12 <- flux_from_nc("MIROC-ESM-CHEM", "2p0", "8p5")
flux_2p0K_13 <- flux_from_nc("NorESM1-M",      "2p0", "8p5")
toc()

tic()
prec_2p0K_1  <- prec_from_nc("HadGEM2-ES",     "2p0", "2p6")
prec_2p0K_2  <- prec_from_nc("IPSL-CM5A-LR",   "2p0", "2p6")
prec_2p0K_3  <- prec_from_nc("MIROC-ESM-CHEM", "2p0", "2p6")
prec_2p0K_4  <- prec_from_nc("GFDL-ESM2M",     "2p0", "6p0")
prec_2p0K_5  <- prec_from_nc("HadGEM2-ES",     "2p0", "6p0")
prec_2p0K_6  <- prec_from_nc("IPSL-CM5A-LR",   "2p0", "6p0")
prec_2p0K_7  <- prec_from_nc("MIROC-ESM-CHEM", "2p0", "6p0")
prec_2p0K_8  <- prec_from_nc("NorESM1-M",      "2p0", "6p0")
prec_2p0K_9  <- prec_from_nc("GFDL-ESM2M",     "2p0", "8p5")
prec_2p0K_10 <- prec_from_nc("HadGEM2-ES",     "2p0", "8p5")
prec_2p0K_11 <- prec_from_nc("IPSL-CM5A-LR",   "2p0", "8p5")
prec_2p0K_12 <- prec_from_nc("MIROC-ESM-CHEM", "2p0", "8p5")
prec_2p0K_13 <- prec_from_nc("NorESM1-M",      "2p0", "8p5")
toc()

#3K warming level
tic()
flux_3p0K_1 <- flux_from_nc("HadGEM2-ES",     "3p0", "6p0")
flux_3p0K_2 <- flux_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
flux_3p0K_3 <- flux_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
flux_3p0K_4 <- flux_from_nc("GFDL-ESM2M",     "3p0", "8p5")
flux_3p0K_5 <- flux_from_nc("HadGEM2-ES",     "3p0", "8p5")
flux_3p0K_6 <- flux_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
flux_3p0K_7 <- flux_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
flux_3p0K_8 <- flux_from_nc("NorESM1-M",      "3p0", "8p5")
toc()

tic()
prec_3p0K_1 <- prec_from_nc("HadGEM2-ES",     "3p0", "6p0")
prec_3p0K_2 <- prec_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
prec_3p0K_3 <- prec_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
prec_3p0K_4 <- prec_from_nc("GFDL-ESM2M",     "3p0", "8p5")
prec_3p0K_5 <- prec_from_nc("HadGEM2-ES",     "3p0", "8p5")
prec_3p0K_6 <- prec_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
prec_3p0K_7 <- prec_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
prec_3p0K_8 <- prec_from_nc("NorESM1-M",      "3p0", "8p5")
toc()

#ann_max----

#get annual maxima discharge characteristics
f_max_mag_doy <- function(data_in, break_day = 274){
  
  my_date <- data_in$date
  start_y <- format(my_date[1], "%Y")
  end_y <- format(my_date[length(my_date)], "%Y")
  
  #remove date column
  data_in <- data_in[-which(colnames(data_in) == "date")]
  
  max_out <- NULL
  for(i in 1:ncol(data_in)){
    
    data_day <- meltimr::ord_day(data_in = data_in[, i],
                                 date = my_date,
                                 start_y = start_y,
                                 end_y = end_y,
                                 break_day = break_day)
    
    val_max <- apply(data_day, 1, max_na)
    
    max_tim <- function(data_in){
      
      doy_max <- which(data_in == max_na(data_in))[1]
      
      return(doy_max)
      
    }
    
    doy_max <- apply(data_day, 1, max_tim)
    
    max_out <- cbind(max_out, val_max, doy_max)
    
  }
  
  return(max_out)
  
}

#Annual maxima: Historical
dmax_hist_1 <- f_max_mag_doy(disc_hist_1) ; fmax_hist_1 <- f_max_mag_doy(flux_hist_1)
dmax_hist_2 <- f_max_mag_doy(disc_hist_2) ; fmax_hist_2 <- f_max_mag_doy(flux_hist_2)
dmax_hist_3 <- f_max_mag_doy(disc_hist_3) ; fmax_hist_3 <- f_max_mag_doy(flux_hist_3)
dmax_hist_4 <- f_max_mag_doy(disc_hist_4) ; fmax_hist_4 <- f_max_mag_doy(flux_hist_4)
dmax_hist_5 <- f_max_mag_doy(disc_hist_5) ; fmax_hist_5 <- f_max_mag_doy(flux_hist_5)

pmax_hist_1 <- f_max_mag_doy(prec_hist_1)
pmax_hist_2 <- f_max_mag_doy(prec_hist_2)
pmax_hist_3 <- f_max_mag_doy(prec_hist_3)
pmax_hist_4 <- f_max_mag_doy(prec_hist_4)
pmax_hist_5 <- f_max_mag_doy(prec_hist_5)

#Annual maxima: 1.5K warming
dmax_1p5K_1  <- f_max_mag_doy(disc_1p5K_1)  ; fmax_1p5K_1  <- f_max_mag_doy(flux_1p5K_1)
dmax_1p5K_2  <- f_max_mag_doy(disc_1p5K_2)  ; fmax_1p5K_2  <- f_max_mag_doy(flux_1p5K_2)
dmax_1p5K_3  <- f_max_mag_doy(disc_1p5K_3)  ; fmax_1p5K_3  <- f_max_mag_doy(flux_1p5K_3)
dmax_1p5K_4  <- f_max_mag_doy(disc_1p5K_4)  ; fmax_1p5K_4  <- f_max_mag_doy(flux_1p5K_4)
dmax_1p5K_5  <- f_max_mag_doy(disc_1p5K_5)  ; fmax_1p5K_5  <- f_max_mag_doy(flux_1p5K_5)
dmax_1p5K_6  <- f_max_mag_doy(disc_1p5K_6)  ; fmax_1p5K_6  <- f_max_mag_doy(flux_1p5K_6)
dmax_1p5K_7  <- f_max_mag_doy(disc_1p5K_7)  ; fmax_1p5K_7  <- f_max_mag_doy(flux_1p5K_7)
dmax_1p5K_8  <- f_max_mag_doy(disc_1p5K_8)  ; fmax_1p5K_8  <- f_max_mag_doy(flux_1p5K_8)
dmax_1p5K_9  <- f_max_mag_doy(disc_1p5K_9)  ; fmax_1p5K_9  <- f_max_mag_doy(flux_1p5K_9)
dmax_1p5K_10 <- f_max_mag_doy(disc_1p5K_10) ; fmax_1p5K_10 <- f_max_mag_doy(flux_1p5K_10)
dmax_1p5K_11 <- f_max_mag_doy(disc_1p5K_11) ; fmax_1p5K_11 <- f_max_mag_doy(flux_1p5K_11)
dmax_1p5K_12 <- f_max_mag_doy(disc_1p5K_12) ; fmax_1p5K_12 <- f_max_mag_doy(flux_1p5K_12)
dmax_1p5K_13 <- f_max_mag_doy(disc_1p5K_13) ; fmax_1p5K_13 <- f_max_mag_doy(flux_1p5K_13)
dmax_1p5K_14 <- f_max_mag_doy(disc_1p5K_14) ; fmax_1p5K_14 <- f_max_mag_doy(flux_1p5K_14)

pmax_1p5K_1  <- f_max_mag_doy(prec_1p5K_1)
pmax_1p5K_2  <- f_max_mag_doy(prec_1p5K_2)
pmax_1p5K_3  <- f_max_mag_doy(prec_1p5K_3)
pmax_1p5K_4  <- f_max_mag_doy(prec_1p5K_4)
pmax_1p5K_5  <- f_max_mag_doy(prec_1p5K_5)
pmax_1p5K_6  <- f_max_mag_doy(prec_1p5K_6)
pmax_1p5K_7  <- f_max_mag_doy(prec_1p5K_7)
pmax_1p5K_8  <- f_max_mag_doy(prec_1p5K_8)
pmax_1p5K_9  <- f_max_mag_doy(prec_1p5K_9)
pmax_1p5K_10 <- f_max_mag_doy(prec_1p5K_10)
pmax_1p5K_11 <- f_max_mag_doy(prec_1p5K_11)
pmax_1p5K_12 <- f_max_mag_doy(prec_1p5K_12)
pmax_1p5K_13 <- f_max_mag_doy(prec_1p5K_13)
pmax_1p5K_14 <- f_max_mag_doy(prec_1p5K_14)

#Annual maxima: 2K warming
dmax_2p0K_1 <-  f_max_mag_doy(disc_2p0K_1)  ; fmax_2p0K_1 <-  f_max_mag_doy(flux_2p0K_1)
dmax_2p0K_2 <-  f_max_mag_doy(disc_2p0K_2)  ; fmax_2p0K_2 <-  f_max_mag_doy(flux_2p0K_2)
dmax_2p0K_3 <-  f_max_mag_doy(disc_2p0K_3)  ; fmax_2p0K_3 <-  f_max_mag_doy(flux_2p0K_3)
dmax_2p0K_4 <-  f_max_mag_doy(disc_2p0K_4)  ; fmax_2p0K_4 <-  f_max_mag_doy(flux_2p0K_4)
dmax_2p0K_5 <-  f_max_mag_doy(disc_2p0K_5)  ; fmax_2p0K_5 <-  f_max_mag_doy(flux_2p0K_5)
dmax_2p0K_6 <-  f_max_mag_doy(disc_2p0K_6)  ; fmax_2p0K_6 <-  f_max_mag_doy(flux_2p0K_6)
dmax_2p0K_7 <-  f_max_mag_doy(disc_2p0K_7)  ; fmax_2p0K_7 <-  f_max_mag_doy(flux_2p0K_7)
dmax_2p0K_8 <-  f_max_mag_doy(disc_2p0K_8)  ; fmax_2p0K_8 <-  f_max_mag_doy(flux_2p0K_8)
dmax_2p0K_9 <-  f_max_mag_doy(disc_2p0K_9)  ; fmax_2p0K_9 <-  f_max_mag_doy(flux_2p0K_9)
dmax_2p0K_10 <- f_max_mag_doy(disc_2p0K_10) ; fmax_2p0K_10 <- f_max_mag_doy(flux_2p0K_10)
dmax_2p0K_11 <- f_max_mag_doy(disc_2p0K_11) ; fmax_2p0K_11 <- f_max_mag_doy(flux_2p0K_11)
dmax_2p0K_12 <- f_max_mag_doy(disc_2p0K_12) ; fmax_2p0K_12 <- f_max_mag_doy(flux_2p0K_12)
dmax_2p0K_13 <- f_max_mag_doy(disc_2p0K_13) ; fmax_2p0K_13 <- f_max_mag_doy(flux_2p0K_13)

pmax_2p0K_1  <- f_max_mag_doy(prec_2p0K_1)
pmax_2p0K_2  <- f_max_mag_doy(prec_2p0K_2)
pmax_2p0K_3  <- f_max_mag_doy(prec_2p0K_3)
pmax_2p0K_4  <- f_max_mag_doy(prec_2p0K_4)
pmax_2p0K_5  <- f_max_mag_doy(prec_2p0K_5)
pmax_2p0K_6  <- f_max_mag_doy(prec_2p0K_6)
pmax_2p0K_7  <- f_max_mag_doy(prec_2p0K_7)
pmax_2p0K_8  <- f_max_mag_doy(prec_2p0K_8)
pmax_2p0K_9  <- f_max_mag_doy(prec_2p0K_9)
pmax_2p0K_10 <- f_max_mag_doy(prec_2p0K_10)
pmax_2p0K_11 <- f_max_mag_doy(prec_2p0K_11)
pmax_2p0K_12 <- f_max_mag_doy(prec_2p0K_12)
pmax_2p0K_13 <- f_max_mag_doy(prec_2p0K_13)

#Annual maxima: 3K warming
dmax_3p0K_1 <- f_max_mag_doy(disc_3p0K_1) ; fmax_3p0K_1 <- f_max_mag_doy(flux_3p0K_1)
dmax_3p0K_2 <- f_max_mag_doy(disc_3p0K_2) ; fmax_3p0K_2 <- f_max_mag_doy(flux_3p0K_2)
dmax_3p0K_3 <- f_max_mag_doy(disc_3p0K_3) ; fmax_3p0K_3 <- f_max_mag_doy(flux_3p0K_3)
dmax_3p0K_4 <- f_max_mag_doy(disc_3p0K_4) ; fmax_3p0K_4 <- f_max_mag_doy(flux_3p0K_4)
dmax_3p0K_5 <- f_max_mag_doy(disc_3p0K_5) ; fmax_3p0K_5 <- f_max_mag_doy(flux_3p0K_5)
dmax_3p0K_6 <- f_max_mag_doy(disc_3p0K_6) ; fmax_3p0K_6 <- f_max_mag_doy(flux_3p0K_6)
dmax_3p0K_7 <- f_max_mag_doy(disc_3p0K_7) ; fmax_3p0K_7 <- f_max_mag_doy(flux_3p0K_7)
dmax_3p0K_8 <- f_max_mag_doy(disc_3p0K_8) ; fmax_3p0K_8 <- f_max_mag_doy(flux_3p0K_8)

pmax_3p0K_1  <- f_max_mag_doy(prec_3p0K_1)
pmax_3p0K_2  <- f_max_mag_doy(prec_3p0K_2)
pmax_3p0K_3  <- f_max_mag_doy(prec_3p0K_3)
pmax_3p0K_4  <- f_max_mag_doy(prec_3p0K_4)
pmax_3p0K_5  <- f_max_mag_doy(prec_3p0K_5)
pmax_3p0K_6  <- f_max_mag_doy(prec_3p0K_6)
pmax_3p0K_7  <- f_max_mag_doy(prec_3p0K_7)
pmax_3p0K_8  <- f_max_mag_doy(prec_3p0K_8)

#Melt fraction peaks function
f_frac_max <- function(flux_in, prec_in, dmax_in){
  
  start_y <- as.numeric(format(flux_in$date[1], "%Y"))
  end_y <- as.numeric(format(flux_in$date[length(flux_in$date)], "%Y"))
  
  #sometimes flux values only starting 2.January
  prec_in <- prec_in[which(prec_in$date %in% flux_in$date), ]
  
  #5 day moving sum total precipitation
  coch_pre_tot_ma <- rollapply(data = prec_in$coch_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  base_pre_tot_ma <- rollapply(data = prec_in$base_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  koel_pre_tot_ma <- rollapply(data = prec_in$koel_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  
  frac_coch <- flux_in$coch_melt_ma / (flux_in$coch_melt_ma + coch_pre_tot_ma)
  frac_base <- flux_in$base_melt_ma / (flux_in$base_melt_ma + base_pre_tot_ma)
  frac_koel <- flux_in$koel_melt_ma / (flux_in$koel_melt_ma + koel_pre_tot_ma)
    
  frac_coch_day <- ord_day(frac_coch,
                           flux_in$date,
                           start_y = start_y,
                           end_y = end_y,
                           break_day = 274)
  
  frac_base_day <- ord_day(frac_base,
                           flux_in$date,
                           start_y = start_y,
                           end_y = end_y,
                           break_day = 274)
  
  frac_koel_day <- ord_day(frac_koel,
                           flux_in$date,
                           start_y = start_y,
                           end_y = end_y,
                           break_day = 274)
  
  frac_max_coch <- NULL
  frac_max_base <- NULL
  frac_max_koel <- NULL
  for(i in 1:nrow(frac_coch_day)){
    frac_max_coch <- c(frac_max_coch, as.numeric(frac_coch_day[i, dmax_in[i, 4]]))
    frac_max_base <- c(frac_max_base, as.numeric(frac_base_day[i, dmax_in[i, 6]]))
    frac_max_koel <- c(frac_max_koel, as.numeric(frac_koel_day[i, dmax_in[i, 6]]))
  }
  
  frac_max_out <- data.frame(coch = frac_max_coch,
                             base = frac_max_base,
                             koel = frac_max_koel)
  
  return(frac_max_out)
  
}

#Melt fraction peaks: Historical
max_frac_hist_1 <- f_frac_max(flux_hist_1, prec_hist_1, dmax_hist_1)
max_frac_hist_2 <- f_frac_max(flux_hist_2, prec_hist_2, dmax_hist_2)
max_frac_hist_3 <- f_frac_max(flux_hist_3, prec_hist_3, dmax_hist_3)
max_frac_hist_4 <- f_frac_max(flux_hist_4, prec_hist_4, dmax_hist_4)
max_frac_hist_5 <- f_frac_max(flux_hist_5, prec_hist_5, dmax_hist_5)

#Melt fraction peaks: 1.5K warming
max_frac_1p5K_1  <- f_frac_max(flux_1p5K_1, prec_1p5K_1, dmax_1p5K_1)
max_frac_1p5K_2  <- f_frac_max(flux_1p5K_2, prec_1p5K_2, dmax_1p5K_2)
max_frac_1p5K_3  <- f_frac_max(flux_1p5K_3, prec_1p5K_3, dmax_1p5K_3)
max_frac_1p5K_4  <- f_frac_max(flux_1p5K_4, prec_1p5K_4, dmax_1p5K_4)
max_frac_1p5K_5  <- f_frac_max(flux_1p5K_5, prec_1p5K_5, dmax_1p5K_5)
max_frac_1p5K_6  <- f_frac_max(flux_1p5K_6, prec_1p5K_6, dmax_1p5K_6)
max_frac_1p5K_7  <- f_frac_max(flux_1p5K_7, prec_1p5K_7, dmax_1p5K_7)
max_frac_1p5K_8  <- f_frac_max(flux_1p5K_8, prec_1p5K_8, dmax_1p5K_8)
max_frac_1p5K_9  <- f_frac_max(flux_1p5K_9, prec_1p5K_9, dmax_1p5K_9)
max_frac_1p5K_10 <- f_frac_max(flux_1p5K_10, prec_1p5K_10, dmax_1p5K_10)
max_frac_1p5K_11 <- f_frac_max(flux_1p5K_11, prec_1p5K_11, dmax_1p5K_11)
max_frac_1p5K_12 <- f_frac_max(flux_1p5K_12, prec_1p5K_12, dmax_1p5K_12)
max_frac_1p5K_13 <- f_frac_max(flux_1p5K_13, prec_1p5K_13, dmax_1p5K_13)
max_frac_1p5K_14 <- f_frac_max(flux_1p5K_14, prec_1p5K_14, dmax_1p5K_14)

#Melt fraction peaks: 2K warming
max_frac_2p0K_1  <- f_frac_max(flux_2p0K_1, prec_2p0K_1, dmax_2p0K_1)
max_frac_2p0K_2  <- f_frac_max(flux_2p0K_2, prec_2p0K_2, dmax_2p0K_2)
max_frac_2p0K_3  <- f_frac_max(flux_2p0K_3, prec_2p0K_3, dmax_2p0K_3)
max_frac_2p0K_4  <- f_frac_max(flux_2p0K_4, prec_2p0K_4, dmax_2p0K_4)
max_frac_2p0K_5  <- f_frac_max(flux_2p0K_5, prec_2p0K_5, dmax_2p0K_5)
max_frac_2p0K_6  <- f_frac_max(flux_2p0K_6, prec_2p0K_6, dmax_2p0K_6)
max_frac_2p0K_7  <- f_frac_max(flux_2p0K_7, prec_2p0K_7, dmax_2p0K_7)
max_frac_2p0K_8  <- f_frac_max(flux_2p0K_8, prec_2p0K_8, dmax_2p0K_8)
max_frac_2p0K_9  <- f_frac_max(flux_2p0K_9, prec_2p0K_9, dmax_2p0K_9)
max_frac_2p0K_10 <- f_frac_max(flux_2p0K_10, prec_2p0K_10, dmax_2p0K_10)
max_frac_2p0K_11 <- f_frac_max(flux_2p0K_11, prec_2p0K_11, dmax_2p0K_11)
max_frac_2p0K_12 <- f_frac_max(flux_2p0K_12, prec_2p0K_12, dmax_2p0K_12)
max_frac_2p0K_13 <- f_frac_max(flux_2p0K_13, prec_2p0K_13, dmax_2p0K_13)

#Melt fraction peaks: 3K warming
max_frac_3p0K_1  <- f_frac_max(flux_3p0K_1, prec_3p0K_1, dmax_3p0K_1)
max_frac_3p0K_2  <- f_frac_max(flux_3p0K_2, prec_3p0K_2, dmax_3p0K_2)
max_frac_3p0K_3  <- f_frac_max(flux_3p0K_3, prec_3p0K_3, dmax_3p0K_3)
max_frac_3p0K_4  <- f_frac_max(flux_3p0K_4, prec_3p0K_4, dmax_3p0K_4)
max_frac_3p0K_5  <- f_frac_max(flux_3p0K_5, prec_3p0K_5, dmax_3p0K_5)
max_frac_3p0K_6  <- f_frac_max(flux_3p0K_6, prec_3p0K_6, dmax_3p0K_6)
max_frac_3p0K_7  <- f_frac_max(flux_3p0K_7, prec_3p0K_7, dmax_3p0K_7)
max_frac_3p0K_8  <- f_frac_max(flux_3p0K_8, prec_3p0K_8, dmax_3p0K_8)

#Put together Historical
max_dis_mag_hist_koel <- c(dmax_hist_1[, 1],  dmax_hist_2[, 1],  dmax_hist_3[, 1],
                           dmax_hist_4[, 1],  dmax_hist_5[, 1])
max_dis_doy_hist_koel <- c(dmax_hist_1[, 2],  dmax_hist_2[, 2],  dmax_hist_3[, 2],
                           dmax_hist_4[, 2],  dmax_hist_5[, 2])
max_dis_mag_hist_coch <- c(dmax_hist_1[, 3],  dmax_hist_2[, 3],  dmax_hist_3[, 3],
                           dmax_hist_4[, 3],  dmax_hist_5[, 3])
max_dis_doy_hist_coch <- c(dmax_hist_1[, 4],  dmax_hist_2[, 4],  dmax_hist_3[, 4],
                           dmax_hist_4[, 4],  dmax_hist_5[, 4])
max_dis_mag_hist_base <- c(dmax_hist_1[, 5],  dmax_hist_2[, 5],  dmax_hist_3[, 5],
                           dmax_hist_4[, 5],  dmax_hist_5[, 5])
max_dis_doy_hist_base <- c(dmax_hist_1[, 6],  dmax_hist_2[, 6],  dmax_hist_3[, 6],
                           dmax_hist_4[, 6],  dmax_hist_5[, 6])
max_dis_fra_hist_base <- c(max_frac_hist_1$base,  max_frac_hist_2$base,  max_frac_hist_3$base,
                           max_frac_hist_4$base,  max_frac_hist_5$base)
max_dis_fra_hist_coch <- c(max_frac_hist_1$coch,  max_frac_hist_2$coch,  max_frac_hist_3$coch,
                           max_frac_hist_4$coch,  max_frac_hist_5$coch)
max_dis_fra_hist_koel <- c(max_frac_hist_1$koel,  max_frac_hist_2$koel,  max_frac_hist_3$koel,
                           max_frac_hist_4$koel,  max_frac_hist_5$koel)

max_mel_mag_hist_base <- c(fmax_hist_1[, 1],  fmax_hist_2[, 1],  fmax_hist_3[, 1],
                           fmax_hist_4[, 1],  fmax_hist_5[, 1])
max_mel_doy_hist_base <- c(fmax_hist_1[, 2],  fmax_hist_2[, 2],  fmax_hist_3[, 2],
                           fmax_hist_4[, 2],  fmax_hist_5[, 2])
max_mel_mag_hist_coch <- c(fmax_hist_1[, 3],  fmax_hist_2[, 3],  fmax_hist_3[, 3],
                           fmax_hist_4[, 3],  fmax_hist_5[, 3])
max_mel_doy_hist_coch <- c(fmax_hist_1[, 4],  fmax_hist_2[, 4],  fmax_hist_3[, 4],
                           fmax_hist_4[, 4],  fmax_hist_5[, 4])
max_mel_mag_hist_koel <- c(fmax_hist_1[, 5],  fmax_hist_2[, 5],  fmax_hist_3[, 5],
                           fmax_hist_4[, 5],  fmax_hist_5[, 5])
max_mel_doy_hist_koel <- c(fmax_hist_1[, 6],  fmax_hist_2[, 6],  fmax_hist_3[, 6],
                           fmax_hist_4[, 6],  fmax_hist_5[, 6])

max_prt_mag_hist_base <- c(pmax_hist_1[, 1],  pmax_hist_2[, 1],  pmax_hist_3[, 1],
                           pmax_hist_4[, 1],  pmax_hist_5[, 1])
max_prt_doy_hist_base <- c(pmax_hist_1[, 2],  pmax_hist_2[, 2],  pmax_hist_3[, 2],
                           pmax_hist_4[, 2],  pmax_hist_5[, 2])
max_prt_mag_hist_coch <- c(pmax_hist_1[, 3],  pmax_hist_2[, 3],  pmax_hist_3[, 3],
                           pmax_hist_4[, 3],  pmax_hist_5[, 3])
max_prt_doy_hist_coch <- c(pmax_hist_1[, 4],  pmax_hist_2[, 4],  pmax_hist_3[, 4],
                           pmax_hist_4[, 4],  pmax_hist_5[, 4])
max_prt_mag_hist_koel <- c(pmax_hist_1[, 5],  pmax_hist_2[, 5],  pmax_hist_3[, 5],
                           pmax_hist_4[, 5],  pmax_hist_5[, 5])
max_prt_doy_hist_koel <- c(pmax_hist_1[, 6],  pmax_hist_2[, 6],  pmax_hist_3[, 6],
                           pmax_hist_4[, 6],  pmax_hist_5[, 6])

max_prl_mag_hist_base <- c(pmax_hist_1[, 7],  pmax_hist_2[, 7],  pmax_hist_3[, 7],
                           pmax_hist_4[, 7],  pmax_hist_5[, 7])
max_prl_doy_hist_base <- c(pmax_hist_1[, 8],  pmax_hist_2[, 8],  pmax_hist_3[, 8],
                           pmax_hist_4[, 8],  pmax_hist_5[, 8])
max_prl_mag_hist_coch <- c(pmax_hist_1[, 9],  pmax_hist_2[, 9],  pmax_hist_3[, 9],
                           pmax_hist_4[, 9],  pmax_hist_5[, 9])
max_prl_doy_hist_coch <- c(pmax_hist_1[, 10],  pmax_hist_2[, 10],  pmax_hist_3[, 10],
                           pmax_hist_4[, 10],  pmax_hist_5[, 10])
max_prl_mag_hist_koel <- c(pmax_hist_1[, 11],  pmax_hist_2[, 11],  pmax_hist_3[, 11],
                           pmax_hist_4[, 11],  pmax_hist_5[, 11])
max_prl_doy_hist_koel <- c(pmax_hist_1[, 12],  pmax_hist_2[, 12],  pmax_hist_3[, 12],
                           pmax_hist_4[, 12],  pmax_hist_5[, 12])

#Put together 1.5K warming level
max_dis_mag_1p5K_koel <- c(dmax_1p5K_1[, 1],  dmax_1p5K_2[, 1],  dmax_1p5K_3[, 1],
                           dmax_1p5K_4[, 1],  dmax_1p5K_5[, 1],  dmax_1p5K_6[, 1],
                           dmax_1p5K_7[, 1],  dmax_1p5K_8[, 1],  dmax_1p5K_9[, 1], 
                           dmax_1p5K_10[, 1], dmax_1p5K_11[, 1], dmax_1p5K_12[, 1], 
                           dmax_1p5K_13[, 1], dmax_1p5K_14[, 1])
max_dis_doy_1p5K_koel <- c(dmax_1p5K_1[, 2],  dmax_1p5K_2[, 2],  dmax_1p5K_3[, 2],
                           dmax_1p5K_4[, 2],  dmax_1p5K_5[, 2],  dmax_1p5K_6[, 2],
                           dmax_1p5K_7[, 2],  dmax_1p5K_8[, 2],  dmax_1p5K_9[, 2], 
                           dmax_1p5K_10[, 2], dmax_1p5K_11[, 2], dmax_1p5K_12[, 2], 
                           dmax_1p5K_13[, 2], dmax_1p5K_14[, 2])
max_dis_mag_1p5K_coch <- c(dmax_1p5K_1[, 3],  dmax_1p5K_2[, 3],  dmax_1p5K_3[, 3],
                           dmax_1p5K_4[, 3],  dmax_1p5K_5[, 3],  dmax_1p5K_6[, 3],
                           dmax_1p5K_7[, 3],  dmax_1p5K_8[, 3],  dmax_1p5K_9[, 3], 
                           dmax_1p5K_10[, 3], dmax_1p5K_11[, 3], dmax_1p5K_12[, 3], 
                           dmax_1p5K_13[, 3], dmax_1p5K_14[, 3])
max_dis_doy_1p5K_coch <- c(dmax_1p5K_1[, 4],  dmax_1p5K_2[, 4],  dmax_1p5K_3[, 4],
                           dmax_1p5K_4[, 4],  dmax_1p5K_5[, 4],  dmax_1p5K_6[, 4],
                           dmax_1p5K_7[, 4],  dmax_1p5K_8[, 4],  dmax_1p5K_9[, 4], 
                           dmax_1p5K_10[, 4], dmax_1p5K_11[, 4], dmax_1p5K_12[, 4], 
                           dmax_1p5K_13[, 4], dmax_1p5K_14[, 4])
max_dis_mag_1p5K_base <- c(dmax_1p5K_1[, 5],  dmax_1p5K_2[, 5],  dmax_1p5K_3[, 5],
                           dmax_1p5K_4[, 5],  dmax_1p5K_5[, 5],  dmax_1p5K_6[, 5],
                           dmax_1p5K_7[, 5],  dmax_1p5K_8[, 5],  dmax_1p5K_9[, 5], 
                           dmax_1p5K_10[, 5], dmax_1p5K_11[, 5], dmax_1p5K_12[, 5], 
                           dmax_1p5K_13[, 5], dmax_1p5K_14[, 5])
max_dis_doy_1p5K_base <- c(dmax_1p5K_1[, 6],  dmax_1p5K_2[, 6],  dmax_1p5K_3[, 6],
                           dmax_1p5K_4[, 6],  dmax_1p5K_5[, 6],  dmax_1p5K_6[, 6],
                           dmax_1p5K_7[, 6],  dmax_1p5K_8[, 6],  dmax_1p5K_9[, 6], 
                           dmax_1p5K_10[, 6], dmax_1p5K_11[, 6], dmax_1p5K_12[, 6], 
                           dmax_1p5K_13[, 6], dmax_1p5K_14[, 6])
max_dis_fra_1p5K_base <- c(max_frac_1p5K_1$base,  max_frac_1p5K_2$base,  max_frac_1p5K_3$base,
                           max_frac_1p5K_4$base,  max_frac_1p5K_5$base,  max_frac_1p5K_6$base,
                           max_frac_1p5K_7$base,  max_frac_1p5K_8$base,  max_frac_1p5K_9$base, 
                           max_frac_1p5K_10$base, max_frac_1p5K_11$base, max_frac_1p5K_12$base, 
                           max_frac_1p5K_13$base, max_frac_1p5K_14$base)
max_dis_fra_1p5K_coch <- c(max_frac_1p5K_1$coch,  max_frac_1p5K_2$coch,  max_frac_1p5K_3$coch,
                           max_frac_1p5K_4$coch,  max_frac_1p5K_5$coch,  max_frac_1p5K_6$coch,
                           max_frac_1p5K_7$coch,  max_frac_1p5K_8$coch,  max_frac_1p5K_9$coch, 
                           max_frac_1p5K_10$coch, max_frac_1p5K_11$coch, max_frac_1p5K_12$coch, 
                           max_frac_1p5K_13$coch, max_frac_1p5K_14$coch)
max_dis_fra_1p5K_koel <- c(max_frac_1p5K_1$koel,  max_frac_1p5K_2$koel,  max_frac_1p5K_3$koel,
                           max_frac_1p5K_4$koel,  max_frac_1p5K_5$koel,  max_frac_1p5K_6$koel,
                           max_frac_1p5K_7$koel,  max_frac_1p5K_8$koel,  max_frac_1p5K_9$koel, 
                           max_frac_1p5K_10$koel, max_frac_1p5K_11$koel, max_frac_1p5K_12$koel, 
                           max_frac_1p5K_13$koel, max_frac_1p5K_14$koel)

max_mel_mag_1p5K_base <- c(fmax_1p5K_1[, 1],  fmax_1p5K_2[, 1],  fmax_1p5K_3[, 1],
                           fmax_1p5K_4[, 1],  fmax_1p5K_5[, 1],  fmax_1p5K_6[, 1],
                           fmax_1p5K_7[, 1],  fmax_1p5K_8[, 1],  fmax_1p5K_9[, 1], 
                           fmax_1p5K_10[, 1], fmax_1p5K_11[, 1], fmax_1p5K_12[, 1], 
                           fmax_1p5K_13[, 1], fmax_1p5K_14[, 1])
max_mel_doy_1p5K_base <- c(fmax_1p5K_1[, 2],  fmax_1p5K_2[, 2],  fmax_1p5K_3[, 2],
                           fmax_1p5K_4[, 2],  fmax_1p5K_5[, 2],  fmax_1p5K_6[, 2],
                           fmax_1p5K_7[, 2],  fmax_1p5K_8[, 2],  fmax_1p5K_9[, 2], 
                           fmax_1p5K_10[, 2], fmax_1p5K_11[, 2], fmax_1p5K_12[, 2], 
                           fmax_1p5K_13[, 2], fmax_1p5K_14[, 2])
max_mel_mag_1p5K_coch <- c(fmax_1p5K_1[, 3],  fmax_1p5K_2[, 3],  fmax_1p5K_3[, 3],
                           fmax_1p5K_4[, 3],  fmax_1p5K_5[, 3],  fmax_1p5K_6[, 3],
                           fmax_1p5K_7[, 3],  fmax_1p5K_8[, 3],  fmax_1p5K_9[, 3], 
                           fmax_1p5K_10[, 3], fmax_1p5K_11[, 3], fmax_1p5K_12[, 3], 
                           fmax_1p5K_13[, 3], fmax_1p5K_14[, 3])
max_mel_doy_1p5K_coch <- c(fmax_1p5K_1[, 4],  fmax_1p5K_2[, 4],  fmax_1p5K_3[, 4],
                           fmax_1p5K_4[, 4],  fmax_1p5K_5[, 4],  fmax_1p5K_6[, 4],
                           fmax_1p5K_7[, 4],  fmax_1p5K_8[, 4],  fmax_1p5K_9[, 4], 
                           fmax_1p5K_10[, 4], fmax_1p5K_11[, 4], fmax_1p5K_12[, 4], 
                           fmax_1p5K_13[, 4], fmax_1p5K_14[, 4])
max_mel_mag_1p5K_koel <- c(fmax_1p5K_1[, 5],  fmax_1p5K_2[, 5],  fmax_1p5K_3[, 5],
                           fmax_1p5K_4[, 5],  fmax_1p5K_5[, 5],  fmax_1p5K_6[, 5],
                           fmax_1p5K_7[, 5],  fmax_1p5K_8[, 5],  fmax_1p5K_9[, 5], 
                           fmax_1p5K_10[, 5], fmax_1p5K_11[, 5], fmax_1p5K_12[, 5], 
                           fmax_1p5K_13[, 5], fmax_1p5K_14[, 5])
max_mel_doy_1p5K_koel <- c(fmax_1p5K_1[, 6],  fmax_1p5K_2[, 6],  fmax_1p5K_3[, 6],
                           fmax_1p5K_4[, 6],  fmax_1p5K_5[, 6],  fmax_1p5K_6[, 6],
                           fmax_1p5K_7[, 6],  fmax_1p5K_8[, 6],  fmax_1p5K_9[, 6], 
                           fmax_1p5K_10[, 6], fmax_1p5K_11[, 6], fmax_1p5K_12[, 6], 
                           fmax_1p5K_13[, 6], fmax_1p5K_14[, 6])

max_prt_mag_1p5K_base <- c(pmax_1p5K_1[, 1],  pmax_1p5K_2[, 1],  pmax_1p5K_3[, 1],
                           pmax_1p5K_4[, 1],  pmax_1p5K_5[, 1],  pmax_1p5K_6[, 1],
                           pmax_1p5K_7[, 1],  pmax_1p5K_8[, 1],  pmax_1p5K_9[, 1], 
                           pmax_1p5K_10[, 1], pmax_1p5K_11[, 1], pmax_1p5K_12[, 1], 
                           pmax_1p5K_13[, 1], pmax_1p5K_14[, 1])
max_prt_doy_1p5K_base <- c(pmax_1p5K_1[, 2],  pmax_1p5K_2[, 2],  pmax_1p5K_3[, 2],
                           pmax_1p5K_4[, 2],  pmax_1p5K_5[, 2],  pmax_1p5K_6[, 2],
                           pmax_1p5K_7[, 2],  pmax_1p5K_8[, 2],  pmax_1p5K_9[, 2], 
                           pmax_1p5K_10[, 2], pmax_1p5K_11[, 2], pmax_1p5K_12[, 2], 
                           pmax_1p5K_13[, 2], pmax_1p5K_14[, 2])
max_prt_mag_1p5K_coch <- c(pmax_1p5K_1[, 3],  pmax_1p5K_2[, 3],  pmax_1p5K_3[, 3],
                           pmax_1p5K_4[, 3],  pmax_1p5K_5[, 3],  pmax_1p5K_6[, 3],
                           pmax_1p5K_7[, 3],  pmax_1p5K_8[, 3],  pmax_1p5K_9[, 3], 
                           pmax_1p5K_10[, 3], pmax_1p5K_11[, 3], pmax_1p5K_12[, 3], 
                           pmax_1p5K_13[, 3], pmax_1p5K_14[, 3])
max_prt_doy_1p5K_coch <- c(pmax_1p5K_1[, 4],  pmax_1p5K_2[, 4],  pmax_1p5K_3[, 4],
                           pmax_1p5K_4[, 4],  pmax_1p5K_5[, 4],  pmax_1p5K_6[, 4],
                           pmax_1p5K_7[, 4],  pmax_1p5K_8[, 4],  pmax_1p5K_9[, 4], 
                           pmax_1p5K_10[, 4], pmax_1p5K_11[, 4], pmax_1p5K_12[, 4], 
                           pmax_1p5K_13[, 4], pmax_1p5K_14[, 4])
max_prt_mag_1p5K_koel <- c(pmax_1p5K_1[, 5],  pmax_1p5K_2[, 5],  pmax_1p5K_3[, 5],
                           pmax_1p5K_4[, 5],  pmax_1p5K_5[, 5],  pmax_1p5K_6[, 5],
                           pmax_1p5K_7[, 5],  pmax_1p5K_8[, 5],  pmax_1p5K_9[, 5], 
                           pmax_1p5K_10[, 5], pmax_1p5K_11[, 5], pmax_1p5K_12[, 5], 
                           pmax_1p5K_13[, 5], pmax_1p5K_14[, 5])
max_prt_doy_1p5K_koel <- c(pmax_1p5K_1[, 6],  pmax_1p5K_2[, 6],  pmax_1p5K_3[, 6],
                           pmax_1p5K_4[, 6],  pmax_1p5K_5[, 6],  pmax_1p5K_6[, 6],
                           pmax_1p5K_7[, 6],  pmax_1p5K_8[, 6],  pmax_1p5K_9[, 6], 
                           pmax_1p5K_10[, 6], pmax_1p5K_11[, 6], pmax_1p5K_12[, 6], 
                           pmax_1p5K_13[, 6], pmax_1p5K_14[, 6])

max_prl_mag_1p5K_base <- c(pmax_1p5K_1[, 7],  pmax_1p5K_2[, 7],  pmax_1p5K_3[, 7],
                           pmax_1p5K_4[, 7],  pmax_1p5K_5[, 7],  pmax_1p5K_6[, 7],
                           pmax_1p5K_7[, 7],  pmax_1p5K_8[, 7],  pmax_1p5K_9[, 7], 
                           pmax_1p5K_10[, 7], pmax_1p5K_11[, 7], pmax_1p5K_12[, 7], 
                           pmax_1p5K_13[, 7], pmax_1p5K_14[, 7])
max_prl_doy_1p5K_base <- c(pmax_1p5K_1[, 8],  pmax_1p5K_2[, 8],  pmax_1p5K_3[, 8],
                           pmax_1p5K_4[, 8],  pmax_1p5K_5[, 8],  pmax_1p5K_6[, 8],
                           pmax_1p5K_7[, 8],  pmax_1p5K_8[, 8],  pmax_1p5K_9[, 8], 
                           pmax_1p5K_10[, 8], pmax_1p5K_11[, 8], pmax_1p5K_12[, 8], 
                           pmax_1p5K_13[, 8], pmax_1p5K_14[, 8])
max_prl_mag_1p5K_coch <- c(pmax_1p5K_1[, 9],  pmax_1p5K_2[, 9],  pmax_1p5K_3[, 9],
                           pmax_1p5K_4[, 9],  pmax_1p5K_5[, 9],  pmax_1p5K_6[, 9],
                           pmax_1p5K_7[, 9],  pmax_1p5K_8[, 9],  pmax_1p5K_9[, 9], 
                           pmax_1p5K_10[, 9], pmax_1p5K_11[, 9], pmax_1p5K_12[, 9], 
                           pmax_1p5K_13[, 9], pmax_1p5K_14[, 9])
max_prl_doy_1p5K_coch <- c(pmax_1p5K_1[, 10],  pmax_1p5K_2[, 10],  pmax_1p5K_3[, 10],
                           pmax_1p5K_4[, 10],  pmax_1p5K_5[, 10],  pmax_1p5K_6[, 10],
                           pmax_1p5K_7[, 10],  pmax_1p5K_8[, 10],  pmax_1p5K_9[, 10], 
                           pmax_1p5K_10[, 10], pmax_1p5K_11[, 10], pmax_1p5K_12[, 10], 
                           pmax_1p5K_13[, 10], pmax_1p5K_14[, 10])
max_prl_mag_1p5K_koel <- c(pmax_1p5K_1[, 11],  pmax_1p5K_2[, 11],  pmax_1p5K_3[, 11],
                           pmax_1p5K_4[, 11],  pmax_1p5K_5[, 11],  pmax_1p5K_6[, 11],
                           pmax_1p5K_7[, 11],  pmax_1p5K_8[, 11],  pmax_1p5K_9[, 11], 
                           pmax_1p5K_10[, 11], pmax_1p5K_11[, 11], pmax_1p5K_12[, 11], 
                           pmax_1p5K_13[, 11], pmax_1p5K_14[, 11])
max_prl_doy_1p5K_koel <- c(pmax_1p5K_1[, 12],  pmax_1p5K_2[, 12],  pmax_1p5K_3[, 12],
                           pmax_1p5K_4[, 12],  pmax_1p5K_5[, 12],  pmax_1p5K_6[, 12],
                           pmax_1p5K_7[, 12],  pmax_1p5K_8[, 12],  pmax_1p5K_9[, 12], 
                           pmax_1p5K_10[, 12], pmax_1p5K_11[, 12], pmax_1p5K_12[, 12], 
                           pmax_1p5K_13[, 12], pmax_1p5K_14[, 12])

#Put together 2.0K warming level
max_dis_mag_2p0K_koel <- c(dmax_2p0K_1[, 1],  dmax_2p0K_2[, 1],  dmax_2p0K_3[, 1],
                           dmax_2p0K_4[, 1],  dmax_2p0K_5[, 1],  dmax_2p0K_6[, 1],
                           dmax_2p0K_7[, 1],  dmax_2p0K_8[, 1],  dmax_2p0K_9[, 1], 
                           dmax_2p0K_10[, 1], dmax_2p0K_11[, 1], dmax_2p0K_12[, 1], 
                           dmax_2p0K_13[, 1])
max_dis_doy_2p0K_koel <- c(dmax_2p0K_1[, 2],  dmax_2p0K_2[, 2],  dmax_2p0K_3[, 2],
                           dmax_2p0K_4[, 2],  dmax_2p0K_5[, 2],  dmax_2p0K_6[, 2],
                           dmax_2p0K_7[, 2],  dmax_2p0K_8[, 2],  dmax_2p0K_9[, 2], 
                           dmax_2p0K_10[, 2], dmax_2p0K_11[, 2], dmax_2p0K_12[, 2], 
                           dmax_2p0K_13[, 2])
max_dis_mag_2p0K_coch <- c(dmax_2p0K_1[, 3],  dmax_2p0K_2[, 3],  dmax_2p0K_3[, 3],
                           dmax_2p0K_4[, 3],  dmax_2p0K_5[, 3],  dmax_2p0K_6[, 3],
                           dmax_2p0K_7[, 3],  dmax_2p0K_8[, 3],  dmax_2p0K_9[, 3], 
                           dmax_2p0K_10[, 3], dmax_2p0K_11[, 3], dmax_2p0K_12[, 3], 
                           dmax_2p0K_13[, 3])
max_dis_doy_2p0K_coch <- c(dmax_2p0K_1[, 4],  dmax_2p0K_2[, 4],  dmax_2p0K_3[, 4],
                           dmax_2p0K_4[, 4],  dmax_2p0K_5[, 4],  dmax_2p0K_6[, 4],
                           dmax_2p0K_7[, 4],  dmax_2p0K_8[, 4],  dmax_2p0K_9[, 4], 
                           dmax_2p0K_10[, 4], dmax_2p0K_11[, 4], dmax_2p0K_12[, 4], 
                           dmax_2p0K_13[, 4])
max_dis_mag_2p0K_base <- c(dmax_2p0K_1[, 5],  dmax_2p0K_2[, 5],  dmax_2p0K_3[, 5],
                           dmax_2p0K_4[, 5],  dmax_2p0K_5[, 5],  dmax_2p0K_6[, 5],
                           dmax_2p0K_7[, 5],  dmax_2p0K_8[, 5],  dmax_2p0K_9[, 5], 
                           dmax_2p0K_10[, 5], dmax_2p0K_11[, 5], dmax_2p0K_12[, 5], 
                           dmax_2p0K_13[, 5])
max_dis_doy_2p0K_base <- c(dmax_2p0K_1[, 6],  dmax_2p0K_2[, 6],  dmax_2p0K_3[, 6],
                           dmax_2p0K_4[, 6],  dmax_2p0K_5[, 6],  dmax_2p0K_6[, 6],
                           dmax_2p0K_7[, 6],  dmax_2p0K_8[, 6],  dmax_2p0K_9[, 6], 
                           dmax_2p0K_10[, 6], dmax_2p0K_11[, 6], dmax_2p0K_12[, 6], 
                           dmax_2p0K_13[, 6])
max_dis_fra_2p0K_coch <- c(max_frac_2p0K_1$coch,  max_frac_2p0K_2$coch,  max_frac_2p0K_3$coch,
                           max_frac_2p0K_4$coch,  max_frac_2p0K_5$coch,  max_frac_2p0K_6$coch,
                           max_frac_2p0K_7$coch,  max_frac_2p0K_8$coch,  max_frac_2p0K_9$coch, 
                           max_frac_2p0K_10$coch, max_frac_2p0K_11$coch, max_frac_2p0K_12$coch, 
                           max_frac_2p0K_13$coch)
max_dis_fra_2p0K_base <- c(max_frac_2p0K_1$base,  max_frac_2p0K_2$base,  max_frac_2p0K_3$base,
                           max_frac_2p0K_4$base,  max_frac_2p0K_5$base,  max_frac_2p0K_6$base,
                           max_frac_2p0K_7$base,  max_frac_2p0K_8$base,  max_frac_2p0K_9$base, 
                           max_frac_2p0K_10$base, max_frac_2p0K_11$base, max_frac_2p0K_12$base, 
                           max_frac_2p0K_13$base)
max_dis_fra_2p0K_koel <- c(max_frac_2p0K_1$koel,  max_frac_2p0K_2$koel,  max_frac_2p0K_3$koel,
                           max_frac_2p0K_4$koel,  max_frac_2p0K_5$koel,  max_frac_2p0K_6$koel,
                           max_frac_2p0K_7$koel,  max_frac_2p0K_8$koel,  max_frac_2p0K_9$koel, 
                           max_frac_2p0K_10$koel, max_frac_2p0K_11$koel, max_frac_2p0K_12$koel, 
                           max_frac_2p0K_13$koel)

max_mel_mag_2p0K_base <- c(fmax_2p0K_1[, 1],  fmax_2p0K_2[, 1],  fmax_2p0K_3[, 1],
                           fmax_2p0K_4[, 1],  fmax_2p0K_5[, 1],  fmax_2p0K_6[, 1],
                           fmax_2p0K_7[, 1],  fmax_2p0K_8[, 1],  fmax_2p0K_9[, 1], 
                           fmax_2p0K_10[, 1], fmax_2p0K_11[, 1], fmax_2p0K_12[, 1], 
                           fmax_2p0K_13[, 1])
max_mel_doy_2p0K_base <- c(fmax_2p0K_1[, 2],  fmax_2p0K_2[, 2],  fmax_2p0K_3[, 2],
                           fmax_2p0K_4[, 2],  fmax_2p0K_5[, 2],  fmax_2p0K_6[, 2],
                           fmax_2p0K_7[, 2],  fmax_2p0K_8[, 2],  fmax_2p0K_9[, 2], 
                           fmax_2p0K_10[, 2], fmax_2p0K_11[, 2], fmax_2p0K_12[, 2], 
                           fmax_2p0K_13[, 2])
max_mel_mag_2p0K_coch <- c(fmax_2p0K_1[, 3],  fmax_2p0K_2[, 3],  fmax_2p0K_3[, 3],
                           fmax_2p0K_4[, 3],  fmax_2p0K_5[, 3],  fmax_2p0K_6[, 3],
                           fmax_2p0K_7[, 3],  fmax_2p0K_8[, 3],  fmax_2p0K_9[, 3], 
                           fmax_2p0K_10[, 3], fmax_2p0K_11[, 3], fmax_2p0K_12[, 3], 
                           fmax_2p0K_13[, 3])
max_mel_doy_2p0K_coch <- c(fmax_2p0K_1[, 4],  fmax_2p0K_2[, 4],  fmax_2p0K_3[, 4],
                           fmax_2p0K_4[, 4],  fmax_2p0K_5[, 4],  fmax_2p0K_6[, 4],
                           fmax_2p0K_7[, 4],  fmax_2p0K_8[, 4],  fmax_2p0K_9[, 4], 
                           fmax_2p0K_10[, 4], fmax_2p0K_11[, 4], fmax_2p0K_12[, 4], 
                           fmax_2p0K_13[, 4])
max_mel_mag_2p0K_koel <- c(fmax_2p0K_1[, 5],  fmax_2p0K_2[, 5],  fmax_2p0K_3[, 5],
                           fmax_2p0K_4[, 5],  fmax_2p0K_5[, 5],  fmax_2p0K_6[, 5],
                           fmax_2p0K_7[, 5],  fmax_2p0K_8[, 5],  fmax_2p0K_9[, 5], 
                           fmax_2p0K_10[, 5], fmax_2p0K_11[, 5], fmax_2p0K_12[, 5], 
                           fmax_2p0K_13[, 5])
max_mel_doy_2p0K_koel <- c(fmax_2p0K_1[, 6],  fmax_2p0K_2[, 6],  fmax_2p0K_3[, 6],
                           fmax_2p0K_4[, 6],  fmax_2p0K_5[, 6],  fmax_2p0K_6[, 6],
                           fmax_2p0K_7[, 6],  fmax_2p0K_8[, 6],  fmax_2p0K_9[, 6], 
                           fmax_2p0K_10[, 6], fmax_2p0K_11[, 6], fmax_2p0K_12[, 6], 
                           fmax_2p0K_13[, 6])

max_prt_mag_2p0K_base <- c(pmax_2p0K_1[, 1],  pmax_2p0K_2[, 1],  pmax_2p0K_3[, 1],
                           pmax_2p0K_4[, 1],  pmax_2p0K_5[, 1],  pmax_2p0K_6[, 1],
                           pmax_2p0K_7[, 1],  pmax_2p0K_8[, 1],  pmax_2p0K_9[, 1], 
                           pmax_2p0K_10[, 1], pmax_2p0K_11[, 1], pmax_2p0K_12[, 1], 
                           pmax_2p0K_13[, 1])
max_prt_doy_2p0K_base <- c(pmax_2p0K_1[, 2],  pmax_2p0K_2[, 2],  pmax_2p0K_3[, 2],
                           pmax_2p0K_4[, 2],  pmax_2p0K_5[, 2],  pmax_2p0K_6[, 2],
                           pmax_2p0K_7[, 2],  pmax_2p0K_8[, 2],  pmax_2p0K_9[, 2], 
                           pmax_2p0K_10[, 2], pmax_2p0K_11[, 2], pmax_2p0K_12[, 2], 
                           pmax_2p0K_13[, 2])
max_prt_mag_2p0K_coch <- c(pmax_2p0K_1[, 3],  pmax_2p0K_2[, 3],  pmax_2p0K_3[, 3],
                           pmax_2p0K_4[, 3],  pmax_2p0K_5[, 3],  pmax_2p0K_6[, 3],
                           pmax_2p0K_7[, 3],  pmax_2p0K_8[, 3],  pmax_2p0K_9[, 3], 
                           pmax_2p0K_10[, 3], pmax_2p0K_11[, 3], pmax_2p0K_12[, 3], 
                           pmax_2p0K_13[, 3])
max_prt_doy_2p0K_coch <- c(pmax_2p0K_1[, 4],  pmax_2p0K_2[, 4],  pmax_2p0K_3[, 4],
                           pmax_2p0K_4[, 4],  pmax_2p0K_5[, 4],  pmax_2p0K_6[, 4],
                           pmax_2p0K_7[, 4],  pmax_2p0K_8[, 4],  pmax_2p0K_9[, 4], 
                           pmax_2p0K_10[, 4], pmax_2p0K_11[, 4], pmax_2p0K_12[, 4], 
                           pmax_2p0K_13[, 4])
max_prt_mag_2p0K_koel <- c(pmax_2p0K_1[, 5],  pmax_2p0K_2[, 5],  pmax_2p0K_3[, 5],
                           pmax_2p0K_4[, 5],  pmax_2p0K_5[, 5],  pmax_2p0K_6[, 5],
                           pmax_2p0K_7[, 5],  pmax_2p0K_8[, 5],  pmax_2p0K_9[, 5], 
                           pmax_2p0K_10[, 5], pmax_2p0K_11[, 5], pmax_2p0K_12[, 5], 
                           pmax_2p0K_13[, 5])
max_prt_doy_2p0K_koel <- c(pmax_2p0K_1[, 6],  pmax_2p0K_2[, 6],  pmax_2p0K_3[, 6],
                           pmax_2p0K_4[, 6],  pmax_2p0K_5[, 6],  pmax_2p0K_6[, 6],
                           pmax_2p0K_7[, 6],  pmax_2p0K_8[, 6],  pmax_2p0K_9[, 6], 
                           pmax_2p0K_10[, 6], pmax_2p0K_11[, 6], pmax_2p0K_12[, 6], 
                           pmax_2p0K_13[, 6])

max_prl_mag_2p0K_base <- c(pmax_2p0K_1[, 7],  pmax_2p0K_2[, 7],  pmax_2p0K_3[, 7],
                           pmax_2p0K_4[, 7],  pmax_2p0K_5[, 7],  pmax_2p0K_6[, 7],
                           pmax_2p0K_7[, 7],  pmax_2p0K_8[, 7],  pmax_2p0K_9[, 7], 
                           pmax_2p0K_10[, 7], pmax_2p0K_11[, 7], pmax_2p0K_12[, 7], 
                           pmax_2p0K_13[, 7])
max_prl_doy_2p0K_base <- c(pmax_2p0K_1[, 8],  pmax_2p0K_2[, 8],  pmax_2p0K_3[, 8],
                           pmax_2p0K_4[, 8],  pmax_2p0K_5[, 8],  pmax_2p0K_6[, 8],
                           pmax_2p0K_7[, 8],  pmax_2p0K_8[, 8],  pmax_2p0K_9[, 8], 
                           pmax_2p0K_10[, 8], pmax_2p0K_11[, 8], pmax_2p0K_12[, 8], 
                           pmax_2p0K_13[, 8])
max_prl_mag_2p0K_coch <- c(pmax_2p0K_1[, 9],  pmax_2p0K_2[, 9],  pmax_2p0K_3[, 9],
                           pmax_2p0K_4[, 9],  pmax_2p0K_5[, 9],  pmax_2p0K_6[, 9],
                           pmax_2p0K_7[, 9],  pmax_2p0K_8[, 9],  pmax_2p0K_9[, 9], 
                           pmax_2p0K_10[, 9], pmax_2p0K_11[, 9], pmax_2p0K_12[, 9], 
                           pmax_2p0K_13[, 9])
max_prl_doy_2p0K_coch <- c(pmax_2p0K_1[, 10],  pmax_2p0K_2[, 10],  pmax_2p0K_3[, 10],
                           pmax_2p0K_4[, 10],  pmax_2p0K_5[, 10],  pmax_2p0K_6[, 10],
                           pmax_2p0K_7[, 10],  pmax_2p0K_8[, 10],  pmax_2p0K_9[, 10], 
                           pmax_2p0K_10[, 10], pmax_2p0K_11[, 10], pmax_2p0K_12[, 10], 
                           pmax_2p0K_13[, 10])
max_prl_mag_2p0K_koel <- c(pmax_2p0K_1[, 11],  pmax_2p0K_2[, 11],  pmax_2p0K_3[, 11],
                           pmax_2p0K_4[, 11],  pmax_2p0K_5[, 11],  pmax_2p0K_6[, 11],
                           pmax_2p0K_7[, 11],  pmax_2p0K_8[, 11],  pmax_2p0K_9[, 11], 
                           pmax_2p0K_10[, 11], pmax_2p0K_11[, 11], pmax_2p0K_12[, 11], 
                           pmax_2p0K_13[, 11])
max_prl_doy_2p0K_koel <- c(pmax_2p0K_1[, 12],  pmax_2p0K_2[, 12],  pmax_2p0K_3[, 12],
                           pmax_2p0K_4[, 12],  pmax_2p0K_5[, 12],  pmax_2p0K_6[, 12],
                           pmax_2p0K_7[, 12],  pmax_2p0K_8[, 12],  pmax_2p0K_9[, 12], 
                           pmax_2p0K_10[, 12], pmax_2p0K_11[, 12], pmax_2p0K_12[, 12], 
                           pmax_2p0K_13[, 12])

#Put together 3.0K warming level
max_dis_mag_3p0K_koel <- c(dmax_3p0K_1[, 1],  dmax_3p0K_2[, 1],  dmax_3p0K_3[, 1],
                           dmax_3p0K_4[, 1],  dmax_3p0K_5[, 1],  dmax_3p0K_6[, 1],
                           dmax_3p0K_7[, 1],  dmax_3p0K_8[, 1])
max_dis_doy_3p0K_koel <- c(dmax_3p0K_1[, 2],  dmax_3p0K_2[, 2],  dmax_3p0K_3[, 2],
                           dmax_3p0K_4[, 2],  dmax_3p0K_5[, 2],  dmax_3p0K_6[, 2],
                           dmax_3p0K_7[, 2],  dmax_3p0K_8[, 2])
max_dis_mag_3p0K_coch <- c(dmax_3p0K_1[, 3],  dmax_3p0K_2[, 3],  dmax_3p0K_3[, 3],
                           dmax_3p0K_4[, 3],  dmax_3p0K_5[, 3],  dmax_3p0K_6[, 3],
                           dmax_3p0K_7[, 3],  dmax_3p0K_8[, 3])
max_dis_doy_3p0K_coch <- c(dmax_3p0K_1[, 4],  dmax_3p0K_2[, 4],  dmax_3p0K_3[, 4],
                           dmax_3p0K_4[, 4],  dmax_3p0K_5[, 4],  dmax_3p0K_6[, 4],
                           dmax_3p0K_7[, 4],  dmax_3p0K_8[, 4])
max_dis_mag_3p0K_base <- c(dmax_3p0K_1[, 5],  dmax_3p0K_2[, 5],  dmax_3p0K_3[, 5],
                           dmax_3p0K_4[, 5],  dmax_3p0K_5[, 5],  dmax_3p0K_6[, 5],
                           dmax_3p0K_7[, 5],  dmax_3p0K_8[, 5])
max_dis_doy_3p0K_base <- c(dmax_3p0K_1[, 6],  dmax_3p0K_2[, 6],  dmax_3p0K_3[, 6],
                           dmax_3p0K_4[, 6],  dmax_3p0K_5[, 6],  dmax_3p0K_6[, 6],
                           dmax_3p0K_7[, 6],  dmax_3p0K_8[, 6])
max_dis_fra_3p0K_coch <- c(max_frac_3p0K_1$coch,  max_frac_3p0K_2$coch,  max_frac_3p0K_3$coch,
                           max_frac_3p0K_4$coch,  max_frac_3p0K_5$coch,  max_frac_3p0K_6$coch,
                           max_frac_3p0K_7$coch,  max_frac_3p0K_8$coch)
max_dis_fra_3p0K_base <- c(max_frac_3p0K_1$base,  max_frac_3p0K_2$base,  max_frac_3p0K_3$base,
                           max_frac_3p0K_4$base,  max_frac_3p0K_5$base,  max_frac_3p0K_6$base,
                           max_frac_3p0K_7$base,  max_frac_3p0K_8$base)
max_dis_fra_3p0K_koel <- c(max_frac_3p0K_1$koel,  max_frac_3p0K_2$koel,  max_frac_3p0K_3$koel,
                           max_frac_3p0K_4$koel,  max_frac_3p0K_5$koel,  max_frac_3p0K_6$koel,
                           max_frac_3p0K_7$koel,  max_frac_3p0K_8$koel)

max_mel_mag_3p0K_base <- c(fmax_3p0K_1[, 1],  fmax_3p0K_2[, 1],  fmax_3p0K_3[, 1],
                           fmax_3p0K_4[, 1],  fmax_3p0K_5[, 1],  fmax_3p0K_6[, 1],
                           fmax_3p0K_7[, 1],  fmax_3p0K_8[, 1])
max_mel_doy_3p0K_base <- c(fmax_3p0K_1[, 2],  fmax_3p0K_2[, 2],  fmax_3p0K_3[, 2],
                           fmax_3p0K_4[, 2],  fmax_3p0K_5[, 2],  fmax_3p0K_6[, 2],
                           fmax_3p0K_7[, 2],  fmax_3p0K_8[, 2])
max_mel_mag_3p0K_coch <- c(fmax_3p0K_1[, 3],  fmax_3p0K_2[, 3],  fmax_3p0K_3[, 3],
                           fmax_3p0K_4[, 3],  fmax_3p0K_5[, 3],  fmax_3p0K_6[, 3],
                           fmax_3p0K_7[, 3],  fmax_3p0K_8[, 3])
max_mel_doy_3p0K_coch <- c(fmax_3p0K_1[, 4],  fmax_3p0K_2[, 4],  fmax_3p0K_3[, 4],
                           fmax_3p0K_4[, 4],  fmax_3p0K_5[, 4],  fmax_3p0K_6[, 4],
                           fmax_3p0K_7[, 4],  fmax_3p0K_8[, 4])
max_mel_mag_3p0K_koel <- c(fmax_3p0K_1[, 5],  fmax_3p0K_2[, 5],  fmax_3p0K_3[, 5],
                           fmax_3p0K_4[, 5],  fmax_3p0K_5[, 5],  fmax_3p0K_6[, 5],
                           fmax_3p0K_7[, 5],  fmax_3p0K_8[, 5])
max_mel_doy_3p0K_koel <- c(fmax_3p0K_1[, 6],  fmax_3p0K_2[, 6],  fmax_3p0K_3[, 6],
                           fmax_3p0K_4[, 6],  fmax_3p0K_5[, 6],  fmax_3p0K_6[, 6],
                           fmax_3p0K_7[, 6],  fmax_3p0K_8[, 6])

max_prt_mag_3p0K_base <- c(pmax_3p0K_1[, 1],  pmax_3p0K_2[, 1],  pmax_3p0K_3[, 1],
                           pmax_3p0K_4[, 1],  pmax_3p0K_5[, 1],  pmax_3p0K_6[, 1],
                           pmax_3p0K_7[, 1],  pmax_3p0K_8[, 1])
max_prt_doy_3p0K_base <- c(pmax_3p0K_1[, 2],  pmax_3p0K_2[, 2],  pmax_3p0K_3[, 2],
                           pmax_3p0K_4[, 2],  pmax_3p0K_5[, 2],  pmax_3p0K_6[, 2],
                           pmax_3p0K_7[, 2],  pmax_3p0K_8[, 2])
max_prt_mag_3p0K_coch <- c(pmax_3p0K_1[, 3],  pmax_3p0K_2[, 3],  pmax_3p0K_3[, 3],
                           pmax_3p0K_4[, 3],  pmax_3p0K_5[, 3],  pmax_3p0K_6[, 3],
                           pmax_3p0K_7[, 3],  pmax_3p0K_8[, 3])
max_prt_doy_3p0K_coch <- c(pmax_3p0K_1[, 4],  pmax_3p0K_2[, 4],  pmax_3p0K_3[, 4],
                           pmax_3p0K_4[, 4],  pmax_3p0K_5[, 4],  pmax_3p0K_6[, 4],
                           pmax_3p0K_7[, 4],  pmax_3p0K_8[, 4])
max_prt_mag_3p0K_koel <- c(pmax_3p0K_1[, 5],  pmax_3p0K_2[, 5],  pmax_3p0K_3[, 5],
                           pmax_3p0K_4[, 5],  pmax_3p0K_5[, 5],  pmax_3p0K_6[, 5],
                           pmax_3p0K_7[, 5],  pmax_3p0K_8[, 5])
max_prt_doy_3p0K_koel <- c(pmax_3p0K_1[, 6],  pmax_3p0K_2[, 6],  pmax_3p0K_3[, 6],
                           pmax_3p0K_4[, 6],  pmax_3p0K_5[, 6],  pmax_3p0K_6[, 6],
                           pmax_3p0K_7[, 6],  pmax_3p0K_8[, 6])

max_prl_mag_3p0K_base <- c(pmax_3p0K_1[, 7],  pmax_3p0K_2[, 7],  pmax_3p0K_3[, 7],
                           pmax_3p0K_4[, 7],  pmax_3p0K_5[, 7],  pmax_3p0K_6[, 7],
                           pmax_3p0K_7[, 7],  pmax_3p0K_8[, 7])
max_prl_doy_3p0K_base <- c(pmax_3p0K_1[, 8],  pmax_3p0K_2[, 8],  pmax_3p0K_3[, 8],
                           pmax_3p0K_4[, 8],  pmax_3p0K_5[, 8],  pmax_3p0K_6[, 8],
                           pmax_3p0K_7[, 8],  pmax_3p0K_8[, 8])
max_prl_mag_3p0K_coch <- c(pmax_3p0K_1[, 9],  pmax_3p0K_2[, 9],  pmax_3p0K_3[, 9],
                           pmax_3p0K_4[, 9],  pmax_3p0K_5[, 9],  pmax_3p0K_6[, 9],
                           pmax_3p0K_7[, 9],  pmax_3p0K_8[, 9])
max_prl_doy_3p0K_coch <- c(pmax_3p0K_1[, 10],  pmax_3p0K_2[, 10],  pmax_3p0K_3[, 10],
                           pmax_3p0K_4[, 10],  pmax_3p0K_5[, 10],  pmax_3p0K_6[, 10],
                           pmax_3p0K_7[, 10],  pmax_3p0K_8[, 10])
max_prl_mag_3p0K_koel <- c(pmax_3p0K_1[, 11],  pmax_3p0K_2[, 11],  pmax_3p0K_3[, 11],
                           pmax_3p0K_4[, 11],  pmax_3p0K_5[, 11],  pmax_3p0K_6[, 11],
                           pmax_3p0K_7[, 11],  pmax_3p0K_8[, 11])
max_prl_doy_3p0K_koel <- c(pmax_3p0K_1[, 12],  pmax_3p0K_2[, 12],  pmax_3p0K_3[, 12],
                           pmax_3p0K_4[, 12],  pmax_3p0K_5[, 12],  pmax_3p0K_6[, 12],
                           pmax_3p0K_7[, 12],  pmax_3p0K_8[, 12])

#Plot-function boxplots
plot_box <- function(max_hist, max_1p5K, max_2p0K, max_3p0K, calc_ylims = F, ylims_in = c(0, 1),
                     y_lab = "", do_legend = F, legend_pos = "topleft"){

  if(calc_ylims){
    ylims <- c(min_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)),
               max_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)))
  }else{
    ylims <- ylims_in
  }
  
  col_hist <- "steelblue4"
  col_1p5K <- "grey25"
  col_2p0K <- "orange3"
  col_3p0K <- "darkred"
  
  max_df <- data.frame(max_hist = range(max_hist, na.rm = T),
                       max_1p5K = range(max_1p5K, na.rm = T),
                       max_2p0K = range(max_2p0K, na.rm = T),
                       max_3p0K = range(max_3p0K, na.rm = T))
  
  boxplot(max_df, boxfill = NA, border = NA, axes = F, ylim = ylims)
  axis(2, mgp=c(3, 0.35, 0), tck = -0.015, cex.axis = 2.0)
  mtext(y_lab, side = 2, line = 2.3, cex = 1.3)
  grid(nx = 0, ny = 6, col = "grey55", lwd = 0.5)
  if(do_legend){
    legend(legend_pos, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
           col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.3,
           box.lwd = 0.0, box.col = "black", bg = "white")
  }
  box()
  
  boxplot(max_hist, ylim = ylims, col = col_hist, axes = F, xaxt = "n", add = TRUE, 
          at = 1, boxwex = 1.4, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(max_1p5K, ylim = ylims, col = col_1p5K, axes = F, xaxt = "n", add = TRUE, 
          at = 2, boxwex = 1.4, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(max_2p0K, ylim = ylims, col = col_2p0K, axes = F, xaxt = "n", add = TRUE, 
          at = 3, boxwex = 1.4, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  boxplot(max_3p0K, ylim = ylims, col = col_3p0K, axes = F, xaxt = "n", add = TRUE, 
          at = 4, boxwex = 1.4, whisklwd = 2, staplelwd = 2, whisklty = 1, notch = T,
          outpch = 19)
  
}

#Plot-function histrograms
plot_hist <- function(max_hist, max_1p5K, max_2p0K, max_3p0K, n_breaks = 50, y_lab = ""){
  
  breaks <- seq(min_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)),
                max_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)), 
                length.out = n_breaks)
  
  
  ylims <- c(0, max_na(c(hist(max_hist, breaks = breaks, plot = F)$density, 
                         hist(max_1p5K, breaks = breaks, plot = F)$density,
                         hist(max_2p0K, breaks = breaks, plot = F)$density, 
                         hist(max_3p0K, breaks = breaks, plot = F)$density)))
  
  par(mar =c(0.6, 3.0, 0.6, 0.5))
  
  hist(max_hist, freq = F, col = col_hist, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims, yaxs = "i", xaxs = "i")
  hist(max_1p5K, freq = F, col = col_1p5K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims, yaxs = "i", xaxs = "i")
  hist(max_2p0K, freq = F, col = col_2p0K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims, yaxs = "i", xaxs = "i")
  hist(max_3p0K, freq = F, col = col_3p0K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims, yaxs = "i", xaxs = "i")
  axis(1, mgp=c(3, 0.85, 0), tck = -0.065, cex.axis = 2.0)
  
  cex_header <- 1.3
  par(mar = c(0,0,0,0))
  
  plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
  mtext(y_lab, side = 3, line = -2.5, cex = cex_header, adj = 0.5)
  
}

#Plot annual maxima discharge
pdf(paste0(bas_dir,"res_figs/max_dis_fut_raw.pdf"), width = 16, height = 12)

par(family = "serif")

layout(matrix(c(37, 37, 37, 37,
                38, 1, 2, 3,
                38, 1, 2, 3,
                38, 1, 2, 3,
                38, 1, 2, 3,
                38, 1, 2, 3,
                38, 4, 9, 14,
                38, 5, 10, 15,
                38, 6, 11, 16,
                38, 7, 12, 17,
                38, 8, 13, 18,
                38, 19, 20, 21,
                38, 19, 20, 21,
                38, 19, 20, 21,
                38, 19, 20, 21,
                38, 19, 20, 21,
                38, 22, 27, 32,
                38, 23, 28, 33,
                38, 24, 29, 34,
                38, 25, 30, 35,
                38, 26, 31, 36),
              21, 4, byrow = T), widths=c(0.1, rep(1, 3)), heights=c(0.5, rep(1, 9), 1.5, rep(1, 9), 0.5))
# layout.show(n = 20)

par(mar = c(1.9, 3.0, 0.8, 0.5))

plot_box(max_dis_mag_hist_base, max_dis_mag_1p5K_base, max_dis_mag_2p0K_base, max_dis_mag_3p0K_base,
         y_lab = "", calc_ylims = T, do_legend = T)

plot_box(max_dis_mag_hist_coch, max_dis_mag_1p5K_coch, max_dis_mag_2p0K_coch, max_dis_mag_3p0K_coch,
         y_lab = "", calc_ylims = T, do_legend = F)

plot_box(max_dis_mag_hist_koel, max_dis_mag_1p5K_koel, max_dis_mag_2p0K_koel, max_dis_mag_3p0K_koel,
         y_lab = "", calc_ylims = T, do_legend = F)

par(mar =c(0.6, 3.0, 0.6, 0.5))

plot_hist(max_dis_mag_hist_base, max_dis_mag_1p5K_base, max_dis_mag_2p0K_base, max_dis_mag_3p0K_base,
          y_lab = "")

plot_hist(max_dis_mag_hist_coch, max_dis_mag_1p5K_coch, max_dis_mag_2p0K_coch, max_dis_mag_3p0K_coch,
          y_lab = "")

plot_hist(max_dis_mag_hist_koel, max_dis_mag_1p5K_koel, max_dis_mag_2p0K_koel, max_dis_mag_3p0K_koel,
          y_lab = "")



par(mar = c(1.9, 3.0, 0.8, 0.5))

plot_box(max_dis_doy_hist_base, max_dis_doy_1p5K_base, max_dis_doy_2p0K_base, max_dis_doy_3p0K_base,
         y_lab = "", calc_ylims = T, do_legend = F)

plot_box(max_dis_doy_hist_coch, max_dis_doy_1p5K_coch, max_dis_doy_2p0K_coch, max_dis_doy_3p0K_coch,
         y_lab = "", calc_ylims = T, do_legend = F)

plot_box(max_dis_doy_hist_koel, max_dis_doy_1p5K_koel, max_dis_doy_2p0K_koel, max_dis_doy_3p0K_koel,
         y_lab = "", calc_ylims = T, do_legend = F)

par(mar =c(0.6, 3.0, 0.6, 0.5))

plot_hist(max_dis_doy_hist_base, max_dis_doy_1p5K_base, max_dis_doy_2p0K_base, max_dis_doy_3p0K_base,
          y_lab = "")

plot_hist(max_dis_doy_hist_coch, max_dis_doy_1p5K_coch, max_dis_doy_2p0K_coch, max_dis_doy_3p0K_coch,
          y_lab = "")

plot_hist(max_dis_doy_hist_koel, max_dis_doy_1p5K_koel, max_dis_doy_2p0K_koel, max_dis_doy_3p0K_koel,
          y_lab = "")

#Gauging station
cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.55, cex = cex_header+0.2, adj = 0.191)
mtext("b) Cochem",
      side = 3, line = -2.55, cex = cex_header+0.2, adj = 0.525)
mtext("c) Cologne",
      side = 3, line = -2.55, cex = cex_header+0.2, adj = 0.875)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Discharge magnitude [m³/s]",  side = 2, line = -2.5, cex = cex_header, adj = 0.85, outer = T)
mtext("2. Discharge timing [DOY]",  side = 2, line = -2.5, cex = cex_header, adj = 0.15, outer = T)

dev.off()





#Plot annual maxima additional features
pdf(paste0(bas_dir,"res_figs/max_add_fut.pdf"), width = 10, height = 10)

par(family = "serif")

layout(matrix(c(11, 11, 11,
                12, 1, 2,
                12, 3, 4,
                12, 5, 6,
                12, 7, 8,
                12, 9, 10),
              6, 3, byrow = T), widths=c(0.12, 1, 1), heights=c(0.15, rep(1, 5)))
# layout.show(n = 12)

par(mar = c(1.5, 3.0, 0.8, 0.5))

#Melt fraction
plot_box(max_dis_fra_hist_base, max_dis_fra_1p5K_base, max_dis_fra_2p0K_base, max_dis_fra_3p0K_base,
         y_lab = "[-]", do_legend = T, legend_pos = "topright")

plot_box(max_dis_fra_hist_coch, max_dis_fra_1p5K_coch, max_dis_fra_2p0K_coch, max_dis_fra_3p0K_coch,
         y_lab = "")

#Melt magnitude
plot_box(max_mel_mag_hist_base, max_mel_mag_1p5K_base, max_mel_mag_2p0K_base, max_mel_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_mel_mag_hist_coch, max_mel_mag_1p5K_coch, max_mel_mag_2p0K_coch, max_mel_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

#Melt timing
plot_box(max_mel_doy_hist_base, max_mel_doy_1p5K_base, max_mel_doy_2p0K_base, max_mel_doy_3p0K_base,
         y_lab = "[DOY]", ylims_in = c(1, 365))

plot_box(max_mel_doy_hist_coch, max_mel_doy_1p5K_coch, max_mel_doy_2p0K_coch, max_mel_doy_3p0K_coch,
         y_lab = "", ylims_in = c(1, 365))

#Precipitation total magnitude
plot_box(max_prt_mag_hist_base, max_prt_mag_1p5K_base, max_prt_mag_2p0K_base, max_prt_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_prt_mag_hist_coch, max_prt_mag_1p5K_coch, max_prt_mag_2p0K_coch, max_prt_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

#Precipitation liquid magnitude
plot_box(max_prl_mag_hist_base, max_prl_mag_1p5K_base, max_prl_mag_2p0K_base, max_prl_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_prl_mag_hist_coch, max_prl_mag_1p5K_coch, max_prl_mag_2p0K_coch, max_prl_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

#Gauging station
cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.28)
mtext("b) Cochem",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.81)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Melt fract.",  side = 2, line = -2.2, cex = cex_header, adj = 0.94, outer = T)
mtext("2. Melt mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.70, outer = T)
mtext("3. Melt tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.490, outer = T)
mtext("4. Prec. tot.",  side = 2, line = -2.2, cex = cex_header, adj = 0.27, outer = T)
mtext("5. Prec. liq.",  side = 2, line = -2.2, cex = cex_header, adj = 0.04, outer = T)

dev.off()







pdf(paste0(bas_dir,"res_figs/max_box_fut.pdf"), width = 16, height = 22)

par(family = "serif")
par(mar = c(0.8, 3.0, 0.8, 0.5))

layout(matrix(c(rep(28, 4),
                29, 1, 2, 3,
                29, 4, 5, 6,
                29, 7, 8, 9,
                29, 10, 11, 12,
                29, 13, 14, 15,
                29, 16, 17, 18,
                29, 19, 20, 21,
                29, 22, 23, 24,
                29, 25, 26, 27),
              10, 4, byrow = T), widths=c(0.15, rep(1, 3)), heights=c(0.2, rep(1, 9)))
# layout.show(n = 29)

#Discharge magnitude
plot_box(max_dis_mag_hist_base, max_dis_mag_1p5K_base, max_dis_mag_2p0K_base, max_dis_mag_3p0K_base,
         y_lab = expression(paste("[m"^"3", "s"^"-1","]")), calc_ylims = T, do_legend = T)

plot_box(max_dis_mag_hist_coch, max_dis_mag_1p5K_coch, max_dis_mag_2p0K_coch, max_dis_mag_3p0K_coch,
         y_lab = "", calc_ylims = T, do_legend = T)

plot_box(max_dis_mag_hist_koel, max_dis_mag_1p5K_koel, max_dis_mag_2p0K_koel, max_dis_mag_3p0K_koel,
         y_lab = "", calc_ylims = T, do_legend = T)

#Discharge timing
plot_box(max_dis_doy_hist_base, max_dis_doy_1p5K_base, max_dis_doy_2p0K_base, max_dis_doy_3p0K_base,
         y_lab = "DOY", ylims_in = c(1, 365))

plot_box(max_dis_doy_hist_coch, max_dis_doy_1p5K_coch, max_dis_doy_2p0K_coch, max_dis_doy_3p0K_coch,
         y_lab = "", ylims_in = c(1, 365))

plot_box(max_dis_doy_hist_koel, max_dis_doy_1p5K_koel, max_dis_doy_2p0K_koel, max_dis_doy_3p0K_koel,
         y_lab = "", ylims_in = c(1, 365))

#Melt fraction
plot_box(max_dis_fra_hist_base, max_dis_fra_1p5K_base, max_dis_fra_2p0K_base, max_dis_fra_3p0K_base,
         y_lab = "[-]")

plot_box(max_dis_fra_hist_coch, max_dis_fra_1p5K_coch, max_dis_fra_2p0K_coch, max_dis_fra_3p0K_coch,
         y_lab = "")

plot_box(max_dis_fra_hist_koel, max_dis_fra_1p5K_koel, max_dis_fra_2p0K_koel, max_dis_fra_3p0K_koel,
         y_lab = "")

#Melt magnitude
plot_box(max_mel_mag_hist_base, max_mel_mag_1p5K_base, max_mel_mag_2p0K_base, max_mel_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_mel_mag_hist_coch, max_mel_mag_1p5K_coch, max_mel_mag_2p0K_coch, max_mel_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

plot_box(max_mel_mag_hist_koel, max_mel_mag_1p5K_koel, max_mel_mag_2p0K_koel, max_mel_mag_3p0K_koel,
         y_lab = "", calc_ylims = T)

#Melt timing
plot_box(max_mel_doy_hist_base, max_mel_doy_1p5K_base, max_mel_doy_2p0K_base, max_mel_doy_3p0K_base,
         y_lab = "DOY", ylims_in = c(1, 365))

plot_box(max_mel_doy_hist_coch, max_mel_doy_1p5K_coch, max_mel_doy_2p0K_coch, max_mel_doy_3p0K_coch,
         y_lab = "", ylims_in = c(1, 365))

plot_box(max_mel_doy_hist_koel, max_mel_doy_1p5K_koel, max_mel_doy_2p0K_koel, max_mel_doy_3p0K_koel,
         y_lab = "", ylims_in = c(1, 365))

#Precipitation total magnitude
plot_box(max_prt_mag_hist_base, max_prt_mag_1p5K_base, max_prt_mag_2p0K_base, max_prt_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_prt_mag_hist_coch, max_prt_mag_1p5K_coch, max_prt_mag_2p0K_coch, max_prt_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

plot_box(max_prt_mag_hist_koel, max_prt_mag_1p5K_koel, max_prt_mag_2p0K_koel, max_prt_mag_3p0K_koel,
         y_lab = "", calc_ylims = T)

#Precipitation total timing
plot_box(max_prt_doy_hist_base, max_prt_doy_1p5K_base, max_prt_doy_2p0K_base, max_prt_doy_3p0K_base,
         y_lab = "DOY", ylims_in = c(1, 365))

plot_box(max_prt_doy_hist_coch, max_prt_doy_1p5K_coch, max_prt_doy_2p0K_coch, max_prt_doy_3p0K_coch,
         y_lab = "", ylims_in = c(1, 365))

plot_box(max_prt_doy_hist_koel, max_prt_doy_1p5K_koel, max_prt_doy_2p0K_koel, max_prt_doy_3p0K_koel,
         y_lab = "", ylims_in = c(1, 365))

#Precipitation liquid magnitude
plot_box(max_prl_mag_hist_base, max_prl_mag_1p5K_base, max_prl_mag_2p0K_base, max_prl_mag_3p0K_base,
         y_lab = "[mm]", calc_ylims = T)

plot_box(max_prl_mag_hist_coch, max_prl_mag_1p5K_coch, max_prl_mag_2p0K_coch, max_prl_mag_3p0K_coch,
         y_lab = "", calc_ylims = T)

plot_box(max_prl_mag_hist_koel, max_prl_mag_1p5K_koel, max_prl_mag_2p0K_koel, max_prl_mag_3p0K_koel,
         y_lab = "", calc_ylims = T)

#Precipitation liquid timing
plot_box(max_prl_doy_hist_base, max_prl_doy_1p5K_base, max_prl_doy_2p0K_base, max_prl_doy_3p0K_base,
         y_lab = "DOY", ylims_in = c(1, 365))

plot_box(max_prl_doy_hist_coch, max_prl_doy_1p5K_coch, max_prl_doy_2p0K_coch, max_prl_doy_3p0K_coch,
         y_lab = "", ylims_in = c(1, 365))

plot_box(max_prl_doy_hist_koel, max_prl_doy_1p5K_koel, max_prl_doy_2p0K_koel, max_prl_doy_3p0K_koel,
         y_lab = "", ylims_in = c(1, 365))

med_na(max_prl_mag_hist_coch)

#Gauging station
cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.191)
mtext("b) Cochem",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.525)
mtext("c) Cologne",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.875)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Disc. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.958, outer = T)
mtext("2. Disc. tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.842, outer = T)
mtext("3. Melt fraction",  side = 2, line = -2.2, cex = cex_header, adj = 0.727, outer = T)
mtext("4. Melt mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.607, outer = T)
mtext("5. Melt tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.495, outer = T)
mtext("6. Prec. tot. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.374, outer = T)
mtext("7. Prec. tot. tim",  side = 2, line = -2.2, cex = cex_header, adj = 0.254, outer = T)
mtext("8. Prec. liq. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.134, outer = T)
mtext("9. Prec. liq. tim",  side = 2, line = -2.2, cex = cex_header, adj = 0.018, outer = T)

dev.off()





plot_hist <- function(max_hist, max_1p5K, max_2p0K, max_3p0K, n_breaks = 50, y_lab = ""){
  
  breaks <- seq(min_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)),
                max_na(c(max_hist, max_1p5K, max_2p0K, max_3p0K)), 
                length.out = n_breaks)
  
  ylims <- c(0, max_na(c(hist(max_hist, breaks = breaks, plot = F)$density, 
                         hist(max_1p5K, breaks = breaks, plot = F)$density,
                         hist(max_2p0K, breaks = breaks, plot = F)$density, 
                         hist(max_3p0K, breaks = breaks, plot = F)$density)))
  
  par(mar =c(0.6, 0.5, 0.2, 0.5))
  
  hist(max_hist, freq = F, col = col_hist, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims)
  hist(max_1p5K, freq = F, col = col_1p5K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims)
  hist(max_2p0K, freq = F, col = col_2p0K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims)
  hist(max_3p0K, freq = F, col = col_3p0K, axes = F, breaks = breaks, 
       ylab = "", xlab = "", main = "", ylim = ylims)
  axis(1, mgp=c(3, 0.19, 0), tck = -0.015, cex.axis = 1.4)

  cex_header <- 1.3
  par(mar = c(0,0,0,0))
  
  plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
  mtext(y_lab, side = 3, line = -2.5, cex = cex_header, adj = 0.5)
  
  }

pdf(paste0(bas_dir,"res_figs/max_his_fut.pdf"), width = 16, height = 25)

col_hist <- "steelblue4"
col_1p5K <- "grey25"
col_2p0K <- "orange3"
col_3p0K <- "darkred"
par(family = "serif")

layout(matrix(c(rep(137, 46),
                136, seq(1,     36+9, 1),
                136, seq(37+9,  72+18, 1),
                136, seq(73+18, 108+27, 1)),
              46, 4, byrow = F), widths=c(0.10, rep(1, 3)), heights=c(1, rep(1, 45)))
# layout.show(n = 137)

plot_hist(max_hist = max_dis_mag_hist_base, max_1p5K = max_dis_mag_1p5K_base, 
          max_2p0K = max_dis_mag_2p0K_base, max_3p0K = max_dis_mag_3p0K_base,
          y_lab = expression(paste("[m"^"3", "s"^"-1","]")))

plot_hist(max_dis_doy_hist_base, max_dis_doy_1p5K_base, max_dis_doy_2p0K_base, max_dis_doy_3p0K_base,
          y_lab = "Day of the year")

plot_hist(max_dis_fra_hist_base, max_dis_fra_1p5K_base, max_dis_fra_2p0K_base, max_dis_fra_3p0K_base,
          y_lab = "Fraction [-]")

plot_hist(max_mel_mag_hist_base, max_mel_mag_1p5K_base, max_mel_mag_2p0K_base, max_mel_mag_3p0K_base,
          y_lab = "[mm]")

plot_hist(max_mel_doy_hist_base, max_mel_doy_1p5K_base, max_mel_doy_2p0K_base, max_mel_doy_3p0K_base,
          y_lab = "Day of the year")

plot_hist(max_prt_mag_hist_base, max_prt_mag_1p5K_base, max_prt_mag_2p0K_base, max_prt_mag_3p0K_base,
          y_lab = "[mm]")

plot_hist(max_prt_doy_hist_base, max_prt_doy_1p5K_base, max_prt_doy_2p0K_base, max_prt_doy_3p0K_base,
          y_lab = "Day of the year")

plot_hist(max_prl_mag_hist_base, max_prl_mag_1p5K_base, max_prl_mag_2p0K_base, max_prl_mag_3p0K_base,
          y_lab = "[mm]")

plot_hist(max_prl_doy_hist_base, max_prl_doy_1p5K_base, max_prl_doy_2p0K_base, max_prl_doy_3p0K_base,
          y_lab = "Day of the year")

plot_hist(max_dis_mag_hist_coch, max_dis_mag_1p5K_coch, max_dis_mag_2p0K_coch, max_dis_mag_3p0K_coch,
          y_lab = expression(paste("[m"^"3", "s"^"-1","]")))

plot_hist(max_dis_doy_hist_coch, max_dis_doy_1p5K_coch, max_dis_doy_2p0K_coch, max_dis_doy_3p0K_coch,
          y_lab = "Day of the year")

plot_hist(max_dis_fra_hist_coch, max_dis_fra_1p5K_coch, max_dis_fra_2p0K_coch, max_dis_fra_3p0K_coch,
          y_lab = "Fraction [-]")

plot_hist(max_mel_mag_hist_coch, max_mel_mag_1p5K_coch, max_mel_mag_2p0K_coch, max_mel_mag_3p0K_coch,
          y_lab = "[mm]")

plot_hist(max_mel_doy_hist_coch, max_mel_doy_1p5K_coch, max_mel_doy_2p0K_coch, max_mel_doy_3p0K_coch,
          y_lab = "DOY")

plot_hist(max_prt_mag_hist_coch, max_prt_mag_1p5K_coch, max_prt_mag_2p0K_coch, max_prt_mag_3p0K_coch,
          y_lab = "[mm]")

plot_hist(max_prt_doy_hist_coch, max_prt_doy_1p5K_coch, max_prt_doy_2p0K_coch, max_prt_doy_3p0K_coch,
          y_lab = "Day of the year")

plot_hist(max_prl_mag_hist_coch, max_prl_mag_1p5K_coch, max_prl_mag_2p0K_coch, max_prl_mag_3p0K_coch,
          y_lab = "[mm]")

plot_hist(max_prl_doy_hist_coch, max_prl_doy_1p5K_coch, max_prl_doy_2p0K_coch, max_prl_doy_3p0K_coch,
          y_lab = "Day of the year")

plot_hist(max_dis_mag_hist_koel, max_dis_mag_1p5K_koel, max_dis_mag_2p0K_koel, max_dis_mag_3p0K_koel,
          y_lab = expression(paste("[m"^"3", "s"^"-1","]")))

plot_hist(max_dis_doy_hist_koel, max_dis_doy_1p5K_koel, max_dis_doy_2p0K_koel, max_dis_doy_3p0K_koel,
          y_lab = "Day of the year")

plot_hist(max_dis_fra_hist_koel, max_dis_fra_1p5K_koel, max_dis_fra_2p0K_koel, max_dis_fra_3p0K_koel,
          y_lab = "Fraction [-]")

plot_hist(max_mel_mag_hist_koel, max_mel_mag_1p5K_koel, max_mel_mag_2p0K_koel, max_mel_mag_3p0K_koel,
          y_lab = "[mm]")

plot_hist(max_mel_doy_hist_koel, max_mel_doy_1p5K_koel, max_mel_doy_2p0K_koel, max_mel_doy_3p0K_koel,
          y_lab = "Day of the year")

plot_hist(max_prt_mag_hist_koel, max_prt_mag_1p5K_koel, max_prt_mag_2p0K_koel, max_prt_mag_3p0K_koel,
          y_lab = "[mm]")

plot_hist(max_prt_doy_hist_koel, max_prt_doy_1p5K_koel, max_prt_doy_2p0K_koel, max_prt_doy_3p0K_koel,
          y_lab = "Day of the year")

plot_hist(max_prl_mag_hist_koel, max_prl_mag_1p5K_koel, max_prl_mag_2p0K_koel, max_prl_mag_3p0K_koel,
          y_lab = "[mm]")

plot_hist(max_prl_doy_hist_koel, max_prl_doy_1p5K_koel, max_prl_doy_2p0K_koel, max_prl_doy_3p0K_koel,
          y_lab = "Day of the year")

#Gauging station
cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.161)
mtext("b) Cochem",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.520)
mtext("c) Cologne",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.875)
par(xpd=TRUE)
legend(22, 90, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
       col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.5,
       box.lwd = 0.0, box.col = "white", bg = "white", ncol = 2)
legend(58, 90, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
       col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.5,
       box.lwd = 0.0, box.col = "white", bg = "white", ncol = 2)
legend(94, 90, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
       col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.5,
       box.lwd = 0.0, box.col = "white", bg = "white", ncol = 2)

par(xpd=FALSE)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Disc. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.958, outer = T)
mtext("2. Disc. tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.842, outer = T)
mtext("3. Melt frac.",  side = 2, line = -2.2, cex = cex_header, adj = 0.727, outer = T)
mtext("4. Melt mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.607, outer = T)
mtext("5. Melt tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.490, outer = T)
mtext("6. Prec. tot. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.384, outer = T)
mtext("7. Prec. tot. tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.264, outer = T)
mtext("8. Prec. liq. mag.",  side = 2, line = -2.2, cex = cex_header, adj = 0.154, outer = T)
mtext("9. Prec. liq. tim.",  side = 2, line = -2.2, cex = cex_header, adj = 0.030, outer = T)

dev.off()


#ann_cyc_flux----

#annual cycles melt fraction
annu_flu <- function(flux_in){
  
  date <- flux_in$date
  start_y <- as.numeric(format(date[1], "%Y"))
  end_y <- as.numeric(format(date[length(date)], "%Y"))
  
  mel_base_day <- ord_day(data_in = flux_in$base_melt_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  lpr_base_day <- ord_day(data_in = flux_in$base_lpre_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  fra_base_day <- ord_day(data_in = flux_in$base_frac,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  ele_base_day <- ord_day(data_in = flux_in$base_mel_ele,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  mel_coch_day <- ord_day(data_in = flux_in$coch_melt_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  lpr_coch_day <- ord_day(data_in = flux_in$coch_lpre_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  fra_coch_day <- ord_day(data_in = flux_in$coch_frac,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  ele_coch_day <- ord_day(data_in = flux_in$coch_mel_ele,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  mel_koel_day <- ord_day(data_in = flux_in$koel_melt_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  lpr_koel_day <- ord_day(data_in = flux_in$koel_lpre_ma,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  fra_koel_day <- ord_day(data_in = flux_in$koel_frac,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  ele_koel_day <- ord_day(data_in = flux_in$koel_mel_ele,
                          date = date,
                          start_y = start_y,
                          end_y = end_y,
                          break_day = 274)
  
  mel_base_mea <- apply(mel_base_day, 2, mea_na)
  mel_coch_mea <- apply(mel_coch_day, 2, mea_na)
  mel_koel_mea <- apply(mel_koel_day, 2, mea_na)
  lpr_base_mea <- apply(lpr_base_day, 2, mea_na)
  lpr_coch_mea <- apply(lpr_coch_day, 2, mea_na)
  lpr_koel_mea <- apply(lpr_koel_day, 2, mea_na)
  fra_base_mea <- apply(fra_base_day, 2, mea_na)
  fra_coch_mea <- apply(fra_coch_day, 2, mea_na)
  fra_koel_mea <- apply(fra_koel_day, 2, mea_na)
  ele_base_mea <- apply(ele_base_day, 2, mea_na)
  ele_coch_mea <- apply(ele_coch_day, 2, mea_na)
  ele_koel_mea <- apply(ele_koel_day, 2, mea_na)
  
  fra_mea_out <- data.frame(mel_base_mea = mel_base_mea, 
                            mel_coch_mea = mel_coch_mea,
                            mel_koel_mea = mel_koel_mea,
                            lpr_base_mea = lpr_base_mea,
                            lpr_coch_mea = lpr_coch_mea,
                            lpr_koel_mea = lpr_koel_mea,
                            fra_base_mea = fra_base_mea, 
                            fra_coch_mea = fra_coch_mea,
                            fra_koel_mea = fra_koel_mea,
                            ele_base_mea = ele_base_mea,
                            ele_coch_mea = ele_coch_mea,
                            ele_koel_mea = ele_koel_mea)
  
  return(fra_mea_out)
  
}
annu_pre <- function(prec_in){
  
  date <- prec_in$date
  start_y <- as.numeric(format(date[1], "%Y"))
  end_y <- as.numeric(format(date[length(date)], "%Y"))
  
  base_pre_tot_ms <- rollapply(data = prec_in$base_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  coch_pre_tot_ms <- rollapply(data = prec_in$coch_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  koel_pre_tot_ms <- rollapply(data = prec_in$koel_pre_tot, width = 5,
                               FUN = sum_na, align = "right", fill = NA)
  
  base_pre_tot_day <- ord_day(data_in = base_pre_tot_ms,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  coch_pre_tot_day <- ord_day(data_in = coch_pre_tot_ms,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  koel_pre_tot_day <- ord_day(data_in = koel_pre_tot_ms,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  base_pre_liq_day <- ord_day(data_in = prec_in$base_pre_liq,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  coch_pre_liq_day <- ord_day(data_in = prec_in$coch_pre_liq,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  koel_pre_liq_day <- ord_day(data_in = prec_in$koel_pre_liq,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  base_pro_eff_day <- ord_day(data_in = prec_in$base_pro_eff,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  coch_pro_eff_day <- ord_day(data_in = prec_in$coch_pro_eff,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  koel_pro_eff_day <- ord_day(data_in = prec_in$koel_pro_eff,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  base_pre_tot_mea <- apply(base_pre_tot_day, 2, mea_na)
  coch_pre_tot_mea <- apply(coch_pre_tot_day, 2, mea_na)
  koel_pre_tot_mea <- apply(koel_pre_tot_day, 2, mea_na)
  base_pre_liq_mea <- apply(base_pre_liq_day, 2, mea_na)
  coch_pre_liq_mea <- apply(coch_pre_liq_day, 2, mea_na)
  koel_pre_liq_mea <- apply(koel_pre_liq_day, 2, mea_na)
  base_pro_eff_mea <- apply(base_pro_eff_day, 2, mea_na)
  coch_pro_eff_mea <- apply(coch_pro_eff_day, 2, mea_na)
  koel_pro_eff_mea <- apply(koel_pro_eff_day, 2, mea_na)
  
  fra_mea_out <- data.frame(base_pre_tot = base_pre_tot_mea, 
                            coch_pre_tot = coch_pre_tot_mea,
                            koel_pre_tot = koel_pre_tot_mea,
                            base_pre_liq = base_pre_liq_mea, 
                            coch_pre_liq = coch_pre_liq_mea,
                            koel_pre_liq = koel_pre_liq_mea,
                            base_pro_eff = base_pro_eff_mea, 
                            coch_pro_eff = coch_pro_eff_mea,
                            koel_pro_eff = koel_pro_eff_mea)
  
  return(fra_mea_out)
  
}
annu_dis <- function(disc_in){
  
  date <- disc_in$date
  start_y <- as.numeric(format(date[1], "%Y"))
  end_y <- as.numeric(format(date[length(date)], "%Y"))
  
  base_disc_day <- ord_day(data_in = disc_in$base,
                              date = date,
                              start_y = start_y,
                              end_y = end_y,
                              break_day = 274)
  
  coch_disc_day <- ord_day(data_in = disc_in$coch,
                           date = date,
                           start_y = start_y,
                           end_y = end_y,
                           break_day = 274)
  
  koel_disc_day <- ord_day(data_in = disc_in$koel,
                           date = date,
                           start_y = start_y,
                           end_y = end_y,
                           break_day = 274)
  
  base_dis_mea <- apply(base_disc_day, 2, mea_na)
  coch_dis_mea <- apply(coch_disc_day, 2, mea_na)
  koel_dis_mea <- apply(koel_disc_day, 2, mea_na)
  
  dis_mea_out <- data.frame(base_dis = base_dis_mea, 
                            coch_dis = coch_dis_mea,
                            koel_dis = koel_dis_mea)
  
  return(dis_mea_out)
  
}

#Historical
ann_flu_hist_1 <- annu_flu(flux_hist_1)
ann_flu_hist_2 <- annu_flu(flux_hist_2)
ann_flu_hist_3 <- annu_flu(flux_hist_3)
ann_flu_hist_4 <- annu_flu(flux_hist_4)
ann_flu_hist_5 <- annu_flu(flux_hist_5)

ann_pre_hist_1 <- annu_pre(prec_hist_1)
ann_pre_hist_2 <- annu_pre(prec_hist_2)
ann_pre_hist_3 <- annu_pre(prec_hist_3)
ann_pre_hist_4 <- annu_pre(prec_hist_4)
ann_pre_hist_5 <- annu_pre(prec_hist_5)

ann_dis_hist_1 <- annu_dis(disc_hist_1)
ann_dis_hist_2 <- annu_dis(disc_hist_2)
ann_dis_hist_3 <- annu_dis(disc_hist_3)
ann_dis_hist_4 <- annu_dis(disc_hist_4)
ann_dis_hist_5 <- annu_dis(disc_hist_5)

ann_mel_hist_base_all <- cbind(ann_flu_hist_1$mel_base_mea, ann_flu_hist_2$mel_base_mea, ann_flu_hist_3$mel_base_mea,
                               ann_flu_hist_4$mel_base_mea, ann_flu_hist_5$mel_base_mea)
ann_mel_hist_coch_all <- cbind(ann_flu_hist_1$mel_coch_mea, ann_flu_hist_2$mel_coch_mea, ann_flu_hist_3$mel_coch_mea,
                               ann_flu_hist_4$mel_coch_mea, ann_flu_hist_5$mel_coch_mea)
ann_mel_hist_koel_all <- cbind(ann_flu_hist_1$mel_koel_mea, ann_flu_hist_2$mel_koel_mea, ann_flu_hist_3$mel_koel_mea,
                               ann_flu_hist_4$mel_koel_mea, ann_flu_hist_5$mel_koel_mea)
ann_lpr_hist_base_all <- cbind(ann_flu_hist_1$lpr_base_mea, ann_flu_hist_2$lpr_base_mea, ann_flu_hist_3$lpr_base_mea,
                               ann_flu_hist_4$lpr_base_mea, ann_flu_hist_5$lpr_base_mea)
ann_lpr_hist_coch_all <- cbind(ann_flu_hist_1$lpr_coch_mea, ann_flu_hist_2$lpr_coch_mea, ann_flu_hist_3$lpr_coch_mea,
                               ann_flu_hist_4$lpr_coch_mea, ann_flu_hist_5$lpr_coch_mea)
ann_lpr_hist_koel_all <- cbind(ann_flu_hist_1$lpr_koel_mea, ann_flu_hist_2$lpr_koel_mea, ann_flu_hist_3$lpr_koel_mea,
                               ann_flu_hist_4$lpr_koel_mea, ann_flu_hist_5$lpr_koel_mea)
ann_fra_hist_base_all <- cbind(ann_flu_hist_1$fra_base_mea, ann_flu_hist_2$fra_base_mea, ann_flu_hist_3$fra_base_mea,
                               ann_flu_hist_4$fra_base_mea, ann_flu_hist_5$fra_base_mea)
ann_fra_hist_coch_all <- cbind(ann_flu_hist_1$fra_coch_mea, ann_flu_hist_2$fra_coch_mea, ann_flu_hist_3$fra_coch_mea,
                               ann_flu_hist_4$fra_coch_mea, ann_flu_hist_5$fra_coch_mea)
ann_fra_hist_koel_all <- cbind(ann_flu_hist_1$fra_koel_mea, ann_flu_hist_2$fra_koel_mea, ann_flu_hist_3$fra_koel_mea,
                               ann_flu_hist_4$fra_koel_mea, ann_flu_hist_5$fra_koel_mea)
ann_ele_hist_base_all <- cbind(ann_flu_hist_1$ele_base_mea, ann_flu_hist_2$ele_base_mea, ann_flu_hist_3$ele_base_mea,
                               ann_flu_hist_4$ele_base_mea, ann_flu_hist_5$ele_base_mea)
ann_ele_hist_coch_all <- cbind(ann_flu_hist_1$ele_coch_mea, ann_flu_hist_2$ele_coch_mea, ann_flu_hist_3$ele_coch_mea,
                               ann_flu_hist_4$ele_coch_mea, ann_flu_hist_5$ele_coch_mea)
ann_ele_hist_koel_all <- cbind(ann_flu_hist_1$ele_koel_mea, ann_flu_hist_2$ele_koel_mea, ann_flu_hist_3$ele_koel_mea,
                               ann_flu_hist_4$ele_koel_mea, ann_flu_hist_5$ele_koel_mea)

ann_prt_hist_base_all <- cbind(ann_pre_hist_1$base_pre_tot, ann_pre_hist_2$base_pre_tot, ann_pre_hist_3$base_pre_tot,
                               ann_pre_hist_4$base_pre_tot, ann_pre_hist_5$base_pre_tot)
ann_prt_hist_coch_all <- cbind(ann_pre_hist_1$coch_pre_tot, ann_pre_hist_2$coch_pre_tot, ann_pre_hist_3$coch_pre_tot,
                               ann_pre_hist_4$coch_pre_tot, ann_pre_hist_5$coch_pre_tot)
ann_prt_hist_koel_all <- cbind(ann_pre_hist_1$koel_pre_tot, ann_pre_hist_2$koel_pre_tot, ann_pre_hist_3$koel_pre_tot,
                               ann_pre_hist_4$koel_pre_tot, ann_pre_hist_5$koel_pre_tot)
ann_prl_hist_base_all <- cbind(ann_pre_hist_1$base_pre_liq, ann_pre_hist_2$base_pre_liq, ann_pre_hist_3$base_pre_liq,
                               ann_pre_hist_4$base_pre_liq, ann_pre_hist_5$base_pre_liq)
ann_prl_hist_coch_all <- cbind(ann_pre_hist_1$coch_pre_liq, ann_pre_hist_2$coch_pre_liq, ann_pre_hist_3$coch_pre_liq,
                               ann_pre_hist_4$coch_pre_liq, ann_pre_hist_5$coch_pre_liq)
ann_prl_hist_koel_all <- cbind(ann_pre_hist_1$koel_pre_liq, ann_pre_hist_2$koel_pre_liq, ann_pre_hist_3$koel_pre_liq,
                               ann_pre_hist_4$koel_pre_liq, ann_pre_hist_5$koel_pre_liq)
ann_eff_hist_base_all <- cbind(ann_pre_hist_1$base_pro_eff, ann_pre_hist_2$base_pro_eff, ann_pre_hist_3$base_pro_eff,
                               ann_pre_hist_4$base_pro_eff, ann_pre_hist_5$base_pro_eff)
ann_eff_hist_coch_all <- cbind(ann_pre_hist_1$coch_pro_eff, ann_pre_hist_2$coch_pro_eff, ann_pre_hist_3$coch_pro_eff,
                               ann_pre_hist_4$coch_pro_eff, ann_pre_hist_5$coch_pro_eff)
ann_eff_hist_koel_all <- cbind(ann_pre_hist_1$koel_pro_eff, ann_pre_hist_2$koel_pro_eff, ann_pre_hist_3$koel_pro_eff,
                               ann_pre_hist_4$koel_pro_eff, ann_pre_hist_5$koel_pro_eff)

ann_dis_hist_base_all <- cbind(ann_dis_hist_1$base_dis, ann_dis_hist_2$base_dis, ann_dis_hist_3$base_dis,
                               ann_dis_hist_4$base_dis, ann_dis_hist_5$base_dis)
ann_dis_hist_coch_all <- cbind(ann_dis_hist_1$coch_dis, ann_dis_hist_2$coch_dis, ann_dis_hist_3$coch_dis,
                               ann_dis_hist_4$coch_dis, ann_dis_hist_5$coch_dis)
ann_dis_hist_koel_all <- cbind(ann_dis_hist_1$koel_dis, ann_dis_hist_2$koel_dis, ann_dis_hist_3$koel_dis,
                               ann_dis_hist_4$koel_dis, ann_dis_hist_5$koel_dis)

ann_mel_hist_base <- apply(ann_mel_hist_base_all, 1, mea_na)
ann_mel_hist_coch <- apply(ann_mel_hist_coch_all, 1, mea_na)
ann_mel_hist_koel <- apply(ann_mel_hist_koel_all, 1, mea_na)
ann_lpr_hist_base <- apply(ann_lpr_hist_base_all, 1, mea_na)
ann_lpr_hist_coch <- apply(ann_lpr_hist_coch_all, 1, mea_na)
ann_lpr_hist_koel <- apply(ann_lpr_hist_koel_all, 1, mea_na)
ann_fra_hist_base <- apply(ann_fra_hist_base_all, 1, mea_na)
ann_fra_hist_coch <- apply(ann_fra_hist_coch_all, 1, mea_na)
ann_fra_hist_koel <- apply(ann_fra_hist_koel_all, 1, mea_na)
ann_ele_hist_base <- apply(ann_ele_hist_base_all, 1, mea_na)
ann_ele_hist_coch <- apply(ann_ele_hist_coch_all, 1, mea_na)
ann_ele_hist_koel <- apply(ann_ele_hist_koel_all, 1, mea_na)

ann_prt_hist_base <- apply(ann_prt_hist_base_all, 1, mea_na)
ann_prt_hist_coch <- apply(ann_prt_hist_coch_all, 1, mea_na)
ann_prt_hist_koel <- apply(ann_prt_hist_koel_all, 1, mea_na)
ann_prl_hist_base <- apply(ann_prl_hist_base_all, 1, mea_na)
ann_prl_hist_coch <- apply(ann_prl_hist_coch_all, 1, mea_na)
ann_prl_hist_koel <- apply(ann_prl_hist_koel_all, 1, mea_na)
ann_eff_hist_base <- apply(ann_eff_hist_base_all, 1, mea_na)
ann_eff_hist_coch <- apply(ann_eff_hist_coch_all, 1, mea_na)
ann_eff_hist_koel <- apply(ann_eff_hist_koel_all, 1, mea_na)

ann_dis_hist_base <- apply(ann_dis_hist_base_all, 1, mea_na)
ann_dis_hist_coch <- apply(ann_dis_hist_coch_all, 1, mea_na)
ann_dis_hist_koel <- apply(ann_dis_hist_koel_all, 1, mea_na)

#1.5K warming level
ann_flu_1p5K_1 <-  annu_flu(flux_1p5K_1)
ann_flu_1p5K_2 <-  annu_flu(flux_1p5K_2)
ann_flu_1p5K_3 <-  annu_flu(flux_1p5K_3)
ann_flu_1p5K_4 <-  annu_flu(flux_1p5K_4)
ann_flu_1p5K_5 <-  annu_flu(flux_1p5K_5)
ann_flu_1p5K_6 <-  annu_flu(flux_1p5K_6)
ann_flu_1p5K_7 <-  annu_flu(flux_1p5K_7)
ann_flu_1p5K_8 <-  annu_flu(flux_1p5K_8)
ann_flu_1p5K_9 <-  annu_flu(flux_1p5K_9)
ann_flu_1p5K_10 <- annu_flu(flux_1p5K_10)
ann_flu_1p5K_11 <- annu_flu(flux_1p5K_11)
ann_flu_1p5K_12 <- annu_flu(flux_1p5K_12)
ann_flu_1p5K_13 <- annu_flu(flux_1p5K_13)
ann_flu_1p5K_14 <- annu_flu(flux_1p5K_14)

ann_pre_1p5K_1 <-  annu_pre(prec_1p5K_1)
ann_pre_1p5K_2 <-  annu_pre(prec_1p5K_2)
ann_pre_1p5K_3 <-  annu_pre(prec_1p5K_3)
ann_pre_1p5K_4 <-  annu_pre(prec_1p5K_4)
ann_pre_1p5K_5 <-  annu_pre(prec_1p5K_5)
ann_pre_1p5K_6 <-  annu_pre(prec_1p5K_6)
ann_pre_1p5K_7 <-  annu_pre(prec_1p5K_7)
ann_pre_1p5K_8 <-  annu_pre(prec_1p5K_8)
ann_pre_1p5K_9 <-  annu_pre(prec_1p5K_9)
ann_pre_1p5K_10 <- annu_pre(prec_1p5K_10)
ann_pre_1p5K_11 <- annu_pre(prec_1p5K_11)
ann_pre_1p5K_12 <- annu_pre(prec_1p5K_12)
ann_pre_1p5K_13 <- annu_pre(prec_1p5K_13)
ann_pre_1p5K_14 <- annu_pre(prec_1p5K_14)

ann_dis_1p5K_1 <-  annu_dis(disc_1p5K_1)
ann_dis_1p5K_2 <-  annu_dis(disc_1p5K_2)
ann_dis_1p5K_3 <-  annu_dis(disc_1p5K_3)
ann_dis_1p5K_4 <-  annu_dis(disc_1p5K_4)
ann_dis_1p5K_5 <-  annu_dis(disc_1p5K_5)
ann_dis_1p5K_6 <-  annu_dis(disc_1p5K_6)
ann_dis_1p5K_7 <-  annu_dis(disc_1p5K_7)
ann_dis_1p5K_8 <-  annu_dis(disc_1p5K_8)
ann_dis_1p5K_9 <-  annu_dis(disc_1p5K_9)
ann_dis_1p5K_10 <- annu_dis(disc_1p5K_10)
ann_dis_1p5K_11 <- annu_dis(disc_1p5K_11)
ann_dis_1p5K_12 <- annu_dis(disc_1p5K_12)
ann_dis_1p5K_13 <- annu_dis(disc_1p5K_13)
ann_dis_1p5K_14 <- annu_dis(disc_1p5K_14)

ann_mel_1p5K_base_all <- cbind(ann_flu_1p5K_1$mel_base_mea,  ann_flu_1p5K_2$mel_base_mea,  ann_flu_1p5K_3$mel_base_mea,
                               ann_flu_1p5K_4$mel_base_mea,  ann_flu_1p5K_5$mel_base_mea,  ann_flu_1p5K_6$mel_base_mea,
                               ann_flu_1p5K_7$mel_base_mea,  ann_flu_1p5K_8$mel_base_mea,  ann_flu_1p5K_9$mel_base_mea,
                               ann_flu_1p5K_10$mel_base_mea, ann_flu_1p5K_11$mel_base_mea, ann_flu_1p5K_12$mel_base_mea,
                               ann_flu_1p5K_13$mel_base_mea, ann_flu_1p5K_14$mel_base_mea)
ann_mel_1p5K_coch_all <- cbind(ann_flu_1p5K_1$mel_coch_mea,  ann_flu_1p5K_2$mel_coch_mea,  ann_flu_1p5K_3$mel_coch_mea,
                               ann_flu_1p5K_4$mel_coch_mea,  ann_flu_1p5K_5$mel_coch_mea,  ann_flu_1p5K_6$mel_coch_mea,
                               ann_flu_1p5K_7$mel_coch_mea,  ann_flu_1p5K_8$mel_coch_mea,  ann_flu_1p5K_9$mel_coch_mea,
                               ann_flu_1p5K_10$mel_coch_mea, ann_flu_1p5K_11$mel_coch_mea, ann_flu_1p5K_12$mel_coch_mea,
                               ann_flu_1p5K_13$mel_coch_mea, ann_flu_1p5K_14$mel_coch_mea)
ann_mel_1p5K_koel_all <- cbind(ann_flu_1p5K_1$mel_koel_mea,  ann_flu_1p5K_2$mel_koel_mea,  ann_flu_1p5K_3$mel_koel_mea,
                               ann_flu_1p5K_4$mel_koel_mea,  ann_flu_1p5K_5$mel_koel_mea,  ann_flu_1p5K_6$mel_koel_mea,
                               ann_flu_1p5K_7$mel_koel_mea,  ann_flu_1p5K_8$mel_koel_mea,  ann_flu_1p5K_9$mel_koel_mea,
                               ann_flu_1p5K_10$mel_koel_mea, ann_flu_1p5K_11$mel_koel_mea, ann_flu_1p5K_12$mel_koel_mea,
                               ann_flu_1p5K_13$mel_koel_mea, ann_flu_1p5K_14$mel_koel_mea)
ann_lpr_1p5K_base_all <- cbind(ann_flu_1p5K_1$lpr_base_mea,  ann_flu_1p5K_2$lpr_base_mea,  ann_flu_1p5K_3$lpr_base_mea,
                               ann_flu_1p5K_4$lpr_base_mea,  ann_flu_1p5K_5$lpr_base_mea,  ann_flu_1p5K_6$lpr_base_mea,
                               ann_flu_1p5K_7$lpr_base_mea,  ann_flu_1p5K_8$lpr_base_mea,  ann_flu_1p5K_9$lpr_base_mea,
                               ann_flu_1p5K_10$lpr_base_mea, ann_flu_1p5K_11$lpr_base_mea, ann_flu_1p5K_12$lpr_base_mea,
                               ann_flu_1p5K_13$lpr_base_mea, ann_flu_1p5K_14$lpr_base_mea)
ann_lpr_1p5K_coch_all <- cbind(ann_flu_1p5K_1$lpr_coch_mea,  ann_flu_1p5K_2$lpr_coch_mea,  ann_flu_1p5K_3$lpr_coch_mea,
                               ann_flu_1p5K_4$lpr_coch_mea,  ann_flu_1p5K_5$lpr_coch_mea,  ann_flu_1p5K_6$lpr_coch_mea,
                               ann_flu_1p5K_7$lpr_coch_mea,  ann_flu_1p5K_8$lpr_coch_mea,  ann_flu_1p5K_9$lpr_coch_mea,
                               ann_flu_1p5K_10$lpr_coch_mea, ann_flu_1p5K_11$lpr_coch_mea, ann_flu_1p5K_12$lpr_coch_mea,
                               ann_flu_1p5K_13$lpr_coch_mea, ann_flu_1p5K_14$lpr_coch_mea)
ann_lpr_1p5K_koel_all <- cbind(ann_flu_1p5K_1$lpr_koel_mea,  ann_flu_1p5K_2$lpr_koel_mea,  ann_flu_1p5K_3$lpr_koel_mea,
                               ann_flu_1p5K_4$lpr_koel_mea,  ann_flu_1p5K_5$lpr_koel_mea,  ann_flu_1p5K_6$lpr_koel_mea,
                               ann_flu_1p5K_7$lpr_koel_mea,  ann_flu_1p5K_8$lpr_koel_mea,  ann_flu_1p5K_9$lpr_koel_mea,
                               ann_flu_1p5K_10$lpr_koel_mea, ann_flu_1p5K_11$lpr_koel_mea, ann_flu_1p5K_12$lpr_koel_mea,
                               ann_flu_1p5K_13$lpr_koel_mea, ann_flu_1p5K_14$lpr_koel_mea)
ann_fra_1p5K_base_all <- cbind(ann_flu_1p5K_1$fra_base_mea,  ann_flu_1p5K_2$fra_base_mea,  ann_flu_1p5K_3$fra_base_mea,
                               ann_flu_1p5K_4$fra_base_mea,  ann_flu_1p5K_5$fra_base_mea,  ann_flu_1p5K_6$fra_base_mea,
                               ann_flu_1p5K_7$fra_base_mea,  ann_flu_1p5K_8$fra_base_mea,  ann_flu_1p5K_9$fra_base_mea,
                               ann_flu_1p5K_10$fra_base_mea, ann_flu_1p5K_11$fra_base_mea, ann_flu_1p5K_12$fra_base_mea,
                               ann_flu_1p5K_13$fra_base_mea, ann_flu_1p5K_14$fra_base_mea)
ann_fra_1p5K_coch_all <- cbind(ann_flu_1p5K_1$fra_coch_mea,  ann_flu_1p5K_2$fra_coch_mea,  ann_flu_1p5K_3$fra_coch_mea,
                               ann_flu_1p5K_4$fra_coch_mea,  ann_flu_1p5K_5$fra_coch_mea,  ann_flu_1p5K_6$fra_coch_mea,
                               ann_flu_1p5K_7$fra_coch_mea,  ann_flu_1p5K_8$fra_coch_mea,  ann_flu_1p5K_9$fra_coch_mea,
                               ann_flu_1p5K_10$fra_coch_mea, ann_flu_1p5K_11$fra_coch_mea, ann_flu_1p5K_12$fra_coch_mea,
                               ann_flu_1p5K_13$fra_coch_mea, ann_flu_1p5K_14$fra_coch_mea)
ann_fra_1p5K_koel_all <- cbind(ann_flu_1p5K_1$fra_koel_mea,  ann_flu_1p5K_2$fra_koel_mea,  ann_flu_1p5K_3$fra_koel_mea,
                               ann_flu_1p5K_4$fra_koel_mea,  ann_flu_1p5K_5$fra_koel_mea,  ann_flu_1p5K_6$fra_koel_mea,
                               ann_flu_1p5K_7$fra_koel_mea,  ann_flu_1p5K_8$fra_koel_mea,  ann_flu_1p5K_9$fra_koel_mea,
                               ann_flu_1p5K_10$fra_koel_mea, ann_flu_1p5K_11$fra_koel_mea, ann_flu_1p5K_12$fra_koel_mea,
                               ann_flu_1p5K_13$fra_koel_mea, ann_flu_1p5K_14$fra_koel_mea)
ann_ele_1p5K_base_all <- cbind(ann_flu_1p5K_1$ele_base_mea,  ann_flu_1p5K_2$ele_base_mea,  ann_flu_1p5K_3$ele_base_mea,
                               ann_flu_1p5K_4$ele_base_mea,  ann_flu_1p5K_5$ele_base_mea,  ann_flu_1p5K_6$ele_base_mea,
                               ann_flu_1p5K_7$ele_base_mea,  ann_flu_1p5K_8$ele_base_mea,  ann_flu_1p5K_9$ele_base_mea,
                               ann_flu_1p5K_10$ele_base_mea, ann_flu_1p5K_11$ele_base_mea, ann_flu_1p5K_12$ele_base_mea,
                               ann_flu_1p5K_13$ele_base_mea, ann_flu_1p5K_14$ele_base_mea)
ann_ele_1p5K_coch_all <- cbind(ann_flu_1p5K_1$ele_coch_mea,  ann_flu_1p5K_2$ele_coch_mea,  ann_flu_1p5K_3$ele_coch_mea,
                               ann_flu_1p5K_4$ele_coch_mea,  ann_flu_1p5K_5$ele_coch_mea,  ann_flu_1p5K_6$ele_coch_mea,
                               ann_flu_1p5K_7$ele_coch_mea,  ann_flu_1p5K_8$ele_coch_mea,  ann_flu_1p5K_9$ele_coch_mea,
                               ann_flu_1p5K_10$ele_coch_mea, ann_flu_1p5K_11$ele_coch_mea, ann_flu_1p5K_12$ele_coch_mea,
                               ann_flu_1p5K_13$ele_coch_mea, ann_flu_1p5K_14$ele_coch_mea)
ann_ele_1p5K_koel_all <- cbind(ann_flu_1p5K_1$ele_koel_mea,  ann_flu_1p5K_2$ele_koel_mea,  ann_flu_1p5K_3$ele_koel_mea,
                               ann_flu_1p5K_4$ele_koel_mea,  ann_flu_1p5K_5$ele_koel_mea,  ann_flu_1p5K_6$ele_koel_mea,
                               ann_flu_1p5K_7$ele_koel_mea,  ann_flu_1p5K_8$ele_koel_mea,  ann_flu_1p5K_9$ele_koel_mea,
                               ann_flu_1p5K_10$ele_koel_mea, ann_flu_1p5K_11$ele_koel_mea, ann_flu_1p5K_12$ele_koel_mea,
                               ann_flu_1p5K_13$ele_koel_mea, ann_flu_1p5K_14$ele_koel_mea)

ann_prt_1p5K_base_all <- cbind(ann_pre_1p5K_1$base_pre_tot,  ann_pre_1p5K_2$base_pre_tot,  ann_pre_1p5K_3$base_pre_tot,
                               ann_pre_1p5K_4$base_pre_tot,  ann_pre_1p5K_5$base_pre_tot,  ann_pre_1p5K_6$base_pre_tot,
                               ann_pre_1p5K_7$base_pre_tot,  ann_pre_1p5K_8$base_pre_tot,  ann_pre_1p5K_9$base_pre_tot,
                               ann_pre_1p5K_10$base_pre_tot, ann_pre_1p5K_11$base_pre_tot, ann_pre_1p5K_12$base_pre_tot,
                               ann_pre_1p5K_13$base_pre_tot, ann_pre_1p5K_14$base_pre_tot)
ann_prt_1p5K_coch_all <- cbind(ann_pre_1p5K_1$coch_pre_tot,  ann_pre_1p5K_2$coch_pre_tot,  ann_pre_1p5K_3$coch_pre_tot,
                               ann_pre_1p5K_4$coch_pre_tot,  ann_pre_1p5K_5$coch_pre_tot,  ann_pre_1p5K_6$coch_pre_tot,
                               ann_pre_1p5K_7$coch_pre_tot,  ann_pre_1p5K_8$coch_pre_tot,  ann_pre_1p5K_9$coch_pre_tot,
                               ann_pre_1p5K_10$coch_pre_tot, ann_pre_1p5K_11$coch_pre_tot, ann_pre_1p5K_12$coch_pre_tot,
                               ann_pre_1p5K_13$coch_pre_tot, ann_pre_1p5K_14$coch_pre_tot)
ann_prt_1p5K_koel_all <- cbind(ann_pre_1p5K_1$koel_pre_tot,  ann_pre_1p5K_2$koel_pre_tot,  ann_pre_1p5K_3$koel_pre_tot,
                               ann_pre_1p5K_4$koel_pre_tot,  ann_pre_1p5K_5$koel_pre_tot,  ann_pre_1p5K_6$koel_pre_tot,
                               ann_pre_1p5K_7$koel_pre_tot,  ann_pre_1p5K_8$koel_pre_tot,  ann_pre_1p5K_9$koel_pre_tot,
                               ann_pre_1p5K_10$koel_pre_tot, ann_pre_1p5K_11$koel_pre_tot, ann_pre_1p5K_12$koel_pre_tot,
                               ann_pre_1p5K_13$koel_pre_tot, ann_pre_1p5K_14$koel_pre_tot)
ann_prl_1p5K_base_all <- cbind(ann_pre_1p5K_1$base_pre_liq,  ann_pre_1p5K_2$base_pre_liq,  ann_pre_1p5K_3$base_pre_liq,
                               ann_pre_1p5K_4$base_pre_liq,  ann_pre_1p5K_5$base_pre_liq,  ann_pre_1p5K_6$base_pre_liq,
                               ann_pre_1p5K_7$base_pre_liq,  ann_pre_1p5K_8$base_pre_liq,  ann_pre_1p5K_9$base_pre_liq,
                               ann_pre_1p5K_10$base_pre_liq, ann_pre_1p5K_11$base_pre_liq, ann_pre_1p5K_12$base_pre_liq,
                               ann_pre_1p5K_13$base_pre_liq, ann_pre_1p5K_14$base_pre_liq)
ann_prl_1p5K_coch_all <- cbind(ann_pre_1p5K_1$coch_pre_liq,  ann_pre_1p5K_2$coch_pre_liq,  ann_pre_1p5K_3$coch_pre_liq,
                               ann_pre_1p5K_4$coch_pre_liq,  ann_pre_1p5K_5$coch_pre_liq,  ann_pre_1p5K_6$coch_pre_liq,
                               ann_pre_1p5K_7$coch_pre_liq,  ann_pre_1p5K_8$coch_pre_liq,  ann_pre_1p5K_9$coch_pre_liq,
                               ann_pre_1p5K_10$coch_pre_liq, ann_pre_1p5K_11$coch_pre_liq, ann_pre_1p5K_12$coch_pre_liq,
                               ann_pre_1p5K_13$coch_pre_liq, ann_pre_1p5K_14$coch_pre_liq)
ann_prl_1p5K_koel_all <- cbind(ann_pre_1p5K_1$koel_pre_liq,  ann_pre_1p5K_2$koel_pre_liq,  ann_pre_1p5K_3$koel_pre_liq,
                               ann_pre_1p5K_4$koel_pre_liq,  ann_pre_1p5K_5$koel_pre_liq,  ann_pre_1p5K_6$koel_pre_liq,
                               ann_pre_1p5K_7$koel_pre_liq,  ann_pre_1p5K_8$koel_pre_liq,  ann_pre_1p5K_9$koel_pre_liq,
                               ann_pre_1p5K_10$koel_pre_liq, ann_pre_1p5K_11$koel_pre_liq, ann_pre_1p5K_12$koel_pre_liq,
                               ann_pre_1p5K_13$koel_pre_liq, ann_pre_1p5K_14$koel_pre_liq)
ann_eff_1p5K_base_all <- cbind(ann_pre_1p5K_1$base_pro_eff,  ann_pre_1p5K_2$base_pro_eff,  ann_pre_1p5K_3$base_pro_eff,
                               ann_pre_1p5K_4$base_pro_eff,  ann_pre_1p5K_5$base_pro_eff,  ann_pre_1p5K_6$base_pro_eff,
                               ann_pre_1p5K_7$base_pro_eff,  ann_pre_1p5K_8$base_pro_eff,  ann_pre_1p5K_9$base_pro_eff,
                               ann_pre_1p5K_10$base_pro_eff, ann_pre_1p5K_11$base_pro_eff, ann_pre_1p5K_12$base_pro_eff,
                               ann_pre_1p5K_13$base_pro_eff, ann_pre_1p5K_14$base_pro_eff)
ann_eff_1p5K_coch_all <- cbind(ann_pre_1p5K_1$coch_pro_eff,  ann_pre_1p5K_2$coch_pro_eff,  ann_pre_1p5K_3$coch_pro_eff,
                               ann_pre_1p5K_4$coch_pro_eff,  ann_pre_1p5K_5$coch_pro_eff,  ann_pre_1p5K_6$coch_pro_eff,
                               ann_pre_1p5K_7$coch_pro_eff,  ann_pre_1p5K_8$coch_pro_eff,  ann_pre_1p5K_9$coch_pro_eff,
                               ann_pre_1p5K_10$coch_pro_eff, ann_pre_1p5K_11$coch_pro_eff, ann_pre_1p5K_12$coch_pro_eff,
                               ann_pre_1p5K_13$coch_pro_eff, ann_pre_1p5K_14$coch_pro_eff)
ann_eff_1p5K_koel_all <- cbind(ann_pre_1p5K_1$koel_pro_eff,  ann_pre_1p5K_2$koel_pro_eff,  ann_pre_1p5K_3$koel_pro_eff,
                               ann_pre_1p5K_4$koel_pro_eff,  ann_pre_1p5K_5$koel_pro_eff,  ann_pre_1p5K_6$koel_pro_eff,
                               ann_pre_1p5K_7$koel_pro_eff,  ann_pre_1p5K_8$koel_pro_eff,  ann_pre_1p5K_9$koel_pro_eff,
                               ann_pre_1p5K_10$koel_pro_eff, ann_pre_1p5K_11$koel_pro_eff, ann_pre_1p5K_12$koel_pro_eff,
                               ann_pre_1p5K_13$koel_pro_eff, ann_pre_1p5K_14$koel_pro_eff)

ann_dis_1p5K_base_all <- cbind(ann_dis_1p5K_1$base_dis,  ann_dis_1p5K_2$base_dis,  ann_dis_1p5K_3$base_dis,
                               ann_dis_1p5K_4$base_dis,  ann_dis_1p5K_5$base_dis,  ann_dis_1p5K_6$base_dis,
                               ann_dis_1p5K_7$base_dis,  ann_dis_1p5K_8$base_dis,  ann_dis_1p5K_9$base_dis,
                               ann_dis_1p5K_10$base_dis, ann_dis_1p5K_11$base_dis, ann_dis_1p5K_12$base_dis,
                               ann_dis_1p5K_13$base_dis, ann_dis_1p5K_14$base_dis)
ann_dis_1p5K_coch_all <- cbind(ann_dis_1p5K_1$coch_dis,  ann_dis_1p5K_2$coch_dis,  ann_dis_1p5K_3$coch_dis,
                               ann_dis_1p5K_4$coch_dis,  ann_dis_1p5K_5$coch_dis,  ann_dis_1p5K_6$coch_dis,
                               ann_dis_1p5K_7$coch_dis,  ann_dis_1p5K_8$coch_dis,  ann_dis_1p5K_9$coch_dis,
                               ann_dis_1p5K_10$coch_dis, ann_dis_1p5K_11$coch_dis, ann_dis_1p5K_12$coch_dis,
                               ann_dis_1p5K_13$coch_dis, ann_dis_1p5K_14$coch_dis)
ann_dis_1p5K_koel_all <- cbind(ann_dis_1p5K_1$koel_dis,  ann_dis_1p5K_2$koel_dis,  ann_dis_1p5K_3$koel_dis,
                               ann_dis_1p5K_4$koel_dis,  ann_dis_1p5K_5$koel_dis,  ann_dis_1p5K_6$koel_dis,
                               ann_dis_1p5K_7$koel_dis,  ann_dis_1p5K_8$koel_dis,  ann_dis_1p5K_9$koel_dis,
                               ann_dis_1p5K_10$koel_dis, ann_dis_1p5K_11$koel_dis, ann_dis_1p5K_12$koel_dis,
                               ann_dis_1p5K_13$koel_dis, ann_dis_1p5K_14$koel_dis)

ann_mel_1p5K_base <- apply(ann_mel_1p5K_base_all, 1, mea_na)
ann_mel_1p5K_coch <- apply(ann_mel_1p5K_coch_all, 1, mea_na)
ann_mel_1p5K_koel <- apply(ann_mel_1p5K_koel_all, 1, mea_na)
ann_lpr_1p5K_base <- apply(ann_lpr_1p5K_base_all, 1, mea_na)
ann_lpr_1p5K_coch <- apply(ann_lpr_1p5K_coch_all, 1, mea_na)
ann_lpr_1p5K_koel <- apply(ann_lpr_1p5K_koel_all, 1, mea_na)
ann_fra_1p5K_base <- apply(ann_fra_1p5K_base_all, 1, mea_na)
ann_fra_1p5K_coch <- apply(ann_fra_1p5K_coch_all, 1, mea_na)
ann_fra_1p5K_koel <- apply(ann_fra_1p5K_koel_all, 1, mea_na)
ann_ele_1p5K_base <- apply(ann_ele_1p5K_base_all, 1, mea_na)
ann_ele_1p5K_coch <- apply(ann_ele_1p5K_coch_all, 1, mea_na)
ann_ele_1p5K_koel <- apply(ann_ele_1p5K_koel_all, 1, mea_na)

ann_prt_1p5K_base <- apply(ann_prt_1p5K_base_all, 1, mea_na)
ann_prt_1p5K_coch <- apply(ann_prt_1p5K_coch_all, 1, mea_na)
ann_prt_1p5K_koel <- apply(ann_prt_1p5K_koel_all, 1, mea_na)
ann_prl_1p5K_base <- apply(ann_prl_1p5K_base_all, 1, mea_na)
ann_prl_1p5K_coch <- apply(ann_prl_1p5K_coch_all, 1, mea_na)
ann_prl_1p5K_koel <- apply(ann_prl_1p5K_koel_all, 1, mea_na)
ann_eff_1p5K_base <- apply(ann_eff_1p5K_base_all, 1, mea_na)
ann_eff_1p5K_coch <- apply(ann_eff_1p5K_coch_all, 1, mea_na)
ann_eff_1p5K_koel <- apply(ann_eff_1p5K_koel_all, 1, mea_na)

ann_dis_1p5K_base <- apply(ann_dis_1p5K_base_all, 1, mea_na)
ann_dis_1p5K_coch <- apply(ann_dis_1p5K_coch_all, 1, mea_na)
ann_dis_1p5K_koel <- apply(ann_dis_1p5K_koel_all, 1, mea_na)

#2.0K warming level
ann_flu_2p0K_1 <-  annu_flu(flux_2p0K_1)
ann_flu_2p0K_2 <-  annu_flu(flux_2p0K_2)
ann_flu_2p0K_3 <-  annu_flu(flux_2p0K_3)
ann_flu_2p0K_4 <-  annu_flu(flux_2p0K_4)
ann_flu_2p0K_5 <-  annu_flu(flux_2p0K_5)
ann_flu_2p0K_6 <-  annu_flu(flux_2p0K_6)
ann_flu_2p0K_7 <-  annu_flu(flux_2p0K_7)
ann_flu_2p0K_8 <-  annu_flu(flux_2p0K_8)
ann_flu_2p0K_9 <-  annu_flu(flux_2p0K_9)
ann_flu_2p0K_10 <- annu_flu(flux_2p0K_10)
ann_flu_2p0K_11 <- annu_flu(flux_2p0K_11)
ann_flu_2p0K_12 <- annu_flu(flux_2p0K_12)
ann_flu_2p0K_13 <- annu_flu(flux_2p0K_13)

ann_pre_2p0K_1 <-  annu_pre(prec_2p0K_1)
ann_pre_2p0K_2 <-  annu_pre(prec_2p0K_2)
ann_pre_2p0K_3 <-  annu_pre(prec_2p0K_3)
ann_pre_2p0K_4 <-  annu_pre(prec_2p0K_4)
ann_pre_2p0K_5 <-  annu_pre(prec_2p0K_5)
ann_pre_2p0K_6 <-  annu_pre(prec_2p0K_6)
ann_pre_2p0K_7 <-  annu_pre(prec_2p0K_7)
ann_pre_2p0K_8 <-  annu_pre(prec_2p0K_8)
ann_pre_2p0K_9 <-  annu_pre(prec_2p0K_9)
ann_pre_2p0K_10 <- annu_pre(prec_2p0K_10)
ann_pre_2p0K_11 <- annu_pre(prec_2p0K_11)
ann_pre_2p0K_12 <- annu_pre(prec_2p0K_12)
ann_pre_2p0K_13 <- annu_pre(prec_2p0K_13)

ann_dis_2p0K_1 <-  annu_dis(disc_2p0K_1)
ann_dis_2p0K_2 <-  annu_dis(disc_2p0K_2)
ann_dis_2p0K_3 <-  annu_dis(disc_2p0K_3)
ann_dis_2p0K_4 <-  annu_dis(disc_2p0K_4)
ann_dis_2p0K_5 <-  annu_dis(disc_2p0K_5)
ann_dis_2p0K_6 <-  annu_dis(disc_2p0K_6)
ann_dis_2p0K_7 <-  annu_dis(disc_2p0K_7)
ann_dis_2p0K_8 <-  annu_dis(disc_2p0K_8)
ann_dis_2p0K_9 <-  annu_dis(disc_2p0K_9)
ann_dis_2p0K_10 <- annu_dis(disc_2p0K_10)
ann_dis_2p0K_11 <- annu_dis(disc_2p0K_11)
ann_dis_2p0K_12 <- annu_dis(disc_2p0K_12)
ann_dis_2p0K_13 <- annu_dis(disc_2p0K_13)

ann_mel_2p0K_base_all <- cbind(ann_flu_2p0K_1$mel_base_mea,  ann_flu_2p0K_2$mel_base_mea,  ann_flu_2p0K_3$mel_base_mea,
                               ann_flu_2p0K_4$mel_base_mea,  ann_flu_2p0K_5$mel_base_mea,  ann_flu_2p0K_6$mel_base_mea,
                               ann_flu_2p0K_7$mel_base_mea,  ann_flu_2p0K_8$mel_base_mea,  ann_flu_2p0K_9$mel_base_mea,
                               ann_flu_2p0K_10$mel_base_mea, ann_flu_2p0K_11$mel_base_mea, ann_flu_2p0K_12$mel_base_mea,
                               ann_flu_2p0K_13$mel_base_mea)
ann_mel_2p0K_coch_all <- cbind(ann_flu_2p0K_1$mel_coch_mea,  ann_flu_2p0K_2$mel_coch_mea,  ann_flu_2p0K_3$mel_coch_mea,
                               ann_flu_2p0K_4$mel_coch_mea,  ann_flu_2p0K_5$mel_coch_mea,  ann_flu_2p0K_6$mel_coch_mea,
                               ann_flu_2p0K_7$mel_coch_mea,  ann_flu_2p0K_8$mel_coch_mea,  ann_flu_2p0K_9$mel_coch_mea,
                               ann_flu_2p0K_10$mel_coch_mea, ann_flu_2p0K_11$mel_coch_mea, ann_flu_2p0K_12$mel_coch_mea,
                               ann_flu_2p0K_13$mel_coch_mea)
ann_mel_2p0K_koel_all <- cbind(ann_flu_2p0K_1$mel_koel_mea,  ann_flu_2p0K_2$mel_koel_mea,  ann_flu_2p0K_3$mel_koel_mea,
                               ann_flu_2p0K_4$mel_koel_mea,  ann_flu_2p0K_5$mel_koel_mea,  ann_flu_2p0K_6$mel_koel_mea,
                               ann_flu_2p0K_7$mel_koel_mea,  ann_flu_2p0K_8$mel_koel_mea,  ann_flu_2p0K_9$mel_koel_mea,
                               ann_flu_2p0K_10$mel_koel_mea, ann_flu_2p0K_11$mel_koel_mea, ann_flu_2p0K_12$mel_koel_mea,
                               ann_flu_2p0K_13$mel_koel_mea)
ann_lpr_2p0K_base_all <- cbind(ann_flu_2p0K_1$lpr_base_mea,  ann_flu_2p0K_2$lpr_base_mea,  ann_flu_2p0K_3$lpr_base_mea,
                               ann_flu_2p0K_4$lpr_base_mea,  ann_flu_2p0K_5$lpr_base_mea,  ann_flu_2p0K_6$lpr_base_mea,
                               ann_flu_2p0K_7$lpr_base_mea,  ann_flu_2p0K_8$lpr_base_mea,  ann_flu_2p0K_9$lpr_base_mea,
                               ann_flu_2p0K_10$lpr_base_mea, ann_flu_2p0K_11$lpr_base_mea, ann_flu_2p0K_12$lpr_base_mea,
                               ann_flu_2p0K_13$lpr_base_mea)
ann_lpr_2p0K_coch_all <- cbind(ann_flu_2p0K_1$lpr_coch_mea,  ann_flu_2p0K_2$lpr_coch_mea,  ann_flu_2p0K_3$lpr_coch_mea,
                               ann_flu_2p0K_4$lpr_coch_mea,  ann_flu_2p0K_5$lpr_coch_mea,  ann_flu_2p0K_6$lpr_coch_mea,
                               ann_flu_2p0K_7$lpr_coch_mea,  ann_flu_2p0K_8$lpr_coch_mea,  ann_flu_2p0K_9$lpr_coch_mea,
                               ann_flu_2p0K_10$lpr_coch_mea, ann_flu_2p0K_11$lpr_coch_mea, ann_flu_2p0K_12$lpr_coch_mea,
                               ann_flu_2p0K_13$lpr_coch_mea)
ann_lpr_2p0K_koel_all <- cbind(ann_flu_2p0K_1$lpr_koel_mea,  ann_flu_2p0K_2$lpr_koel_mea,  ann_flu_2p0K_3$lpr_koel_mea,
                               ann_flu_2p0K_4$lpr_koel_mea,  ann_flu_2p0K_5$lpr_koel_mea,  ann_flu_2p0K_6$lpr_koel_mea,
                               ann_flu_2p0K_7$lpr_koel_mea,  ann_flu_2p0K_8$lpr_koel_mea,  ann_flu_2p0K_9$lpr_koel_mea,
                               ann_flu_2p0K_10$lpr_koel_mea, ann_flu_2p0K_11$lpr_koel_mea, ann_flu_2p0K_12$lpr_koel_mea,
                               ann_flu_2p0K_13$lpr_koel_mea)
ann_fra_2p0K_base_all <- cbind(ann_flu_2p0K_1$fra_base_mea,  ann_flu_2p0K_2$fra_base_mea,  ann_flu_2p0K_3$fra_base_mea,
                               ann_flu_2p0K_4$fra_base_mea,  ann_flu_2p0K_5$fra_base_mea,  ann_flu_2p0K_6$fra_base_mea,
                               ann_flu_2p0K_7$fra_base_mea,  ann_flu_2p0K_8$fra_base_mea,  ann_flu_2p0K_9$fra_base_mea,
                               ann_flu_2p0K_10$fra_base_mea, ann_flu_2p0K_11$fra_base_mea, ann_flu_2p0K_12$fra_base_mea,
                               ann_flu_2p0K_13$fra_base_mea)
ann_fra_2p0K_coch_all <- cbind(ann_flu_2p0K_1$fra_coch_mea,  ann_flu_2p0K_2$fra_coch_mea,  ann_flu_2p0K_3$fra_coch_mea,
                               ann_flu_2p0K_4$fra_coch_mea,  ann_flu_2p0K_5$fra_coch_mea,  ann_flu_2p0K_6$fra_coch_mea,
                               ann_flu_2p0K_7$fra_coch_mea,  ann_flu_2p0K_8$fra_coch_mea,  ann_flu_2p0K_9$fra_coch_mea,
                               ann_flu_2p0K_10$fra_coch_mea, ann_flu_2p0K_11$fra_coch_mea, ann_flu_2p0K_12$fra_coch_mea,
                               ann_flu_2p0K_13$fra_coch_mea)
ann_fra_2p0K_koel_all <- cbind(ann_flu_2p0K_1$fra_koel_mea,  ann_flu_2p0K_2$fra_koel_mea,  ann_flu_2p0K_3$fra_koel_mea,
                               ann_flu_2p0K_4$fra_koel_mea,  ann_flu_2p0K_5$fra_koel_mea,  ann_flu_2p0K_6$fra_koel_mea,
                               ann_flu_2p0K_7$fra_koel_mea,  ann_flu_2p0K_8$fra_koel_mea,  ann_flu_2p0K_9$fra_koel_mea,
                               ann_flu_2p0K_10$fra_koel_mea, ann_flu_2p0K_11$fra_koel_mea, ann_flu_2p0K_12$fra_koel_mea,
                               ann_flu_2p0K_13$fra_koel_mea)
ann_ele_2p0K_base_all <- cbind(ann_flu_2p0K_1$ele_base_mea,  ann_flu_2p0K_2$ele_base_mea,  ann_flu_2p0K_3$ele_base_mea,
                               ann_flu_2p0K_4$ele_base_mea,  ann_flu_2p0K_5$ele_base_mea,  ann_flu_2p0K_6$ele_base_mea,
                               ann_flu_2p0K_7$ele_base_mea,  ann_flu_2p0K_8$ele_base_mea,  ann_flu_2p0K_9$ele_base_mea,
                               ann_flu_2p0K_10$ele_base_mea, ann_flu_2p0K_11$ele_base_mea, ann_flu_2p0K_12$ele_base_mea,
                               ann_flu_2p0K_13$ele_base_mea)
ann_ele_2p0K_coch_all <- cbind(ann_flu_2p0K_1$ele_coch_mea,  ann_flu_2p0K_2$ele_coch_mea,  ann_flu_2p0K_3$ele_coch_mea,
                               ann_flu_2p0K_4$ele_coch_mea,  ann_flu_2p0K_5$ele_coch_mea,  ann_flu_2p0K_6$ele_coch_mea,
                               ann_flu_2p0K_7$ele_coch_mea,  ann_flu_2p0K_8$ele_coch_mea,  ann_flu_2p0K_9$ele_coch_mea,
                               ann_flu_2p0K_10$ele_coch_mea, ann_flu_2p0K_11$ele_coch_mea, ann_flu_2p0K_12$ele_coch_mea,
                               ann_flu_2p0K_13$ele_coch_mea)
ann_ele_2p0K_koel_all <- cbind(ann_flu_2p0K_1$ele_koel_mea,  ann_flu_2p0K_2$ele_koel_mea,  ann_flu_2p0K_3$ele_koel_mea,
                               ann_flu_2p0K_4$ele_koel_mea,  ann_flu_2p0K_5$ele_koel_mea,  ann_flu_2p0K_6$ele_koel_mea,
                               ann_flu_2p0K_7$ele_koel_mea,  ann_flu_2p0K_8$ele_koel_mea,  ann_flu_2p0K_9$ele_koel_mea,
                               ann_flu_2p0K_10$ele_koel_mea, ann_flu_2p0K_11$ele_koel_mea, ann_flu_2p0K_12$ele_koel_mea,
                               ann_flu_2p0K_13$ele_koel_mea)

ann_prt_2p0K_base_all <- cbind(ann_pre_2p0K_1$base_pre_tot,  ann_pre_2p0K_2$base_pre_tot,  ann_pre_2p0K_3$base_pre_tot,
                               ann_pre_2p0K_4$base_pre_tot,  ann_pre_2p0K_5$base_pre_tot,  ann_pre_2p0K_6$base_pre_tot,
                               ann_pre_2p0K_7$base_pre_tot,  ann_pre_2p0K_8$base_pre_tot,  ann_pre_2p0K_9$base_pre_tot,
                               ann_pre_2p0K_10$base_pre_tot, ann_pre_2p0K_11$base_pre_tot, ann_pre_2p0K_12$base_pre_tot,
                               ann_pre_2p0K_13$base_pre_tot)
ann_prt_2p0K_coch_all <- cbind(ann_pre_2p0K_1$coch_pre_tot,  ann_pre_2p0K_2$coch_pre_tot,  ann_pre_2p0K_3$coch_pre_tot,
                               ann_pre_2p0K_4$coch_pre_tot,  ann_pre_2p0K_5$coch_pre_tot,  ann_pre_2p0K_6$coch_pre_tot,
                               ann_pre_2p0K_7$coch_pre_tot,  ann_pre_2p0K_8$coch_pre_tot,  ann_pre_2p0K_9$coch_pre_tot,
                               ann_pre_2p0K_10$coch_pre_tot, ann_pre_2p0K_11$coch_pre_tot, ann_pre_2p0K_12$coch_pre_tot,
                               ann_pre_2p0K_13$coch_pre_tot)
ann_prt_2p0K_koel_all <- cbind(ann_pre_2p0K_1$koel_pre_tot,  ann_pre_2p0K_2$koel_pre_tot,  ann_pre_2p0K_3$koel_pre_tot,
                               ann_pre_2p0K_4$koel_pre_tot,  ann_pre_2p0K_5$koel_pre_tot,  ann_pre_2p0K_6$koel_pre_tot,
                               ann_pre_2p0K_7$koel_pre_tot,  ann_pre_2p0K_8$koel_pre_tot,  ann_pre_2p0K_9$koel_pre_tot,
                               ann_pre_2p0K_10$koel_pre_tot, ann_pre_2p0K_11$koel_pre_tot, ann_pre_2p0K_12$koel_pre_tot,
                               ann_pre_2p0K_13$koel_pre_tot)
ann_prl_2p0K_base_all <- cbind(ann_pre_2p0K_1$base_pre_liq,  ann_pre_2p0K_2$base_pre_liq,  ann_pre_2p0K_3$base_pre_liq,
                               ann_pre_2p0K_4$base_pre_liq,  ann_pre_2p0K_5$base_pre_liq,  ann_pre_2p0K_6$base_pre_liq,
                               ann_pre_2p0K_7$base_pre_liq,  ann_pre_2p0K_8$base_pre_liq,  ann_pre_2p0K_9$base_pre_liq,
                               ann_pre_2p0K_10$base_pre_liq, ann_pre_2p0K_11$base_pre_liq, ann_pre_2p0K_12$base_pre_liq,
                               ann_pre_2p0K_13$base_pre_liq)
ann_prl_2p0K_coch_all <- cbind(ann_pre_2p0K_1$coch_pre_liq,  ann_pre_2p0K_2$coch_pre_liq,  ann_pre_2p0K_3$coch_pre_liq,
                               ann_pre_2p0K_4$coch_pre_liq,  ann_pre_2p0K_5$coch_pre_liq,  ann_pre_2p0K_6$coch_pre_liq,
                               ann_pre_2p0K_7$coch_pre_liq,  ann_pre_2p0K_8$coch_pre_liq,  ann_pre_2p0K_9$coch_pre_liq,
                               ann_pre_2p0K_10$coch_pre_liq, ann_pre_2p0K_11$coch_pre_liq, ann_pre_2p0K_12$coch_pre_liq,
                               ann_pre_2p0K_13$coch_pre_liq)
ann_prl_2p0K_koel_all <- cbind(ann_pre_2p0K_1$koel_pre_liq,  ann_pre_2p0K_2$koel_pre_liq,  ann_pre_2p0K_3$koel_pre_liq,
                               ann_pre_2p0K_4$koel_pre_liq,  ann_pre_2p0K_5$koel_pre_liq,  ann_pre_2p0K_6$koel_pre_liq,
                               ann_pre_2p0K_7$koel_pre_liq,  ann_pre_2p0K_8$koel_pre_liq,  ann_pre_2p0K_9$koel_pre_liq,
                               ann_pre_2p0K_10$koel_pre_liq, ann_pre_2p0K_11$koel_pre_liq, ann_pre_2p0K_12$koel_pre_liq,
                               ann_pre_2p0K_13$koel_pre_liq)
ann_eff_2p0K_base_all <- cbind(ann_pre_2p0K_1$base_pro_eff,  ann_pre_2p0K_2$base_pro_eff,  ann_pre_2p0K_3$base_pro_eff,
                               ann_pre_2p0K_4$base_pro_eff,  ann_pre_2p0K_5$base_pro_eff,  ann_pre_2p0K_6$base_pro_eff,
                               ann_pre_2p0K_7$base_pro_eff,  ann_pre_2p0K_8$base_pro_eff,  ann_pre_2p0K_9$base_pro_eff,
                               ann_pre_2p0K_10$base_pro_eff, ann_pre_2p0K_11$base_pro_eff, ann_pre_2p0K_12$base_pro_eff,
                               ann_pre_2p0K_13$base_pro_eff)
ann_eff_2p0K_coch_all <- cbind(ann_pre_2p0K_1$coch_pro_eff,  ann_pre_2p0K_2$coch_pro_eff,  ann_pre_2p0K_3$coch_pro_eff,
                               ann_pre_2p0K_4$coch_pro_eff,  ann_pre_2p0K_5$coch_pro_eff,  ann_pre_2p0K_6$coch_pro_eff,
                               ann_pre_2p0K_7$coch_pro_eff,  ann_pre_2p0K_8$coch_pro_eff,  ann_pre_2p0K_9$coch_pro_eff,
                               ann_pre_2p0K_10$coch_pro_eff, ann_pre_2p0K_11$coch_pro_eff, ann_pre_2p0K_12$coch_pro_eff,
                               ann_pre_2p0K_13$coch_pro_eff)
ann_eff_2p0K_koel_all <- cbind(ann_pre_2p0K_1$koel_pro_eff,  ann_pre_2p0K_2$koel_pro_eff,  ann_pre_2p0K_3$koel_pro_eff,
                               ann_pre_2p0K_4$koel_pro_eff,  ann_pre_2p0K_5$koel_pro_eff,  ann_pre_2p0K_6$koel_pro_eff,
                               ann_pre_2p0K_7$koel_pro_eff,  ann_pre_2p0K_8$koel_pro_eff,  ann_pre_2p0K_9$koel_pro_eff,
                               ann_pre_2p0K_10$koel_pro_eff, ann_pre_2p0K_11$koel_pro_eff, ann_pre_2p0K_12$koel_pro_eff,
                               ann_pre_2p0K_13$koel_pro_eff)

ann_dis_2p0K_base_all <- cbind(ann_dis_2p0K_1$base_dis,  ann_dis_2p0K_2$base_dis,  ann_dis_2p0K_3$base_dis,
                               ann_dis_2p0K_4$base_dis,  ann_dis_2p0K_5$base_dis,  ann_dis_2p0K_6$base_dis,
                               ann_dis_2p0K_7$base_dis,  ann_dis_2p0K_8$base_dis,  ann_dis_2p0K_9$base_dis,
                               ann_dis_2p0K_10$base_dis, ann_dis_2p0K_11$base_dis, ann_dis_2p0K_12$base_dis,
                               ann_dis_2p0K_13$base_dis)
ann_dis_2p0K_coch_all <- cbind(ann_dis_2p0K_1$coch_dis,  ann_dis_2p0K_2$coch_dis,  ann_dis_2p0K_3$coch_dis,
                               ann_dis_2p0K_4$coch_dis,  ann_dis_2p0K_5$coch_dis,  ann_dis_2p0K_6$coch_dis,
                               ann_dis_2p0K_7$coch_dis,  ann_dis_2p0K_8$coch_dis,  ann_dis_2p0K_9$coch_dis,
                               ann_dis_2p0K_10$coch_dis, ann_dis_2p0K_11$coch_dis, ann_dis_2p0K_12$coch_dis,
                               ann_dis_2p0K_13$coch_dis)
ann_dis_2p0K_koel_all <- cbind(ann_dis_2p0K_1$koel_dis,  ann_dis_2p0K_2$koel_dis,  ann_dis_2p0K_3$koel_dis,
                               ann_dis_2p0K_4$koel_dis,  ann_dis_2p0K_5$koel_dis,  ann_dis_2p0K_6$koel_dis,
                               ann_dis_2p0K_7$koel_dis,  ann_dis_2p0K_8$koel_dis,  ann_dis_2p0K_9$koel_dis,
                               ann_dis_2p0K_10$koel_dis, ann_dis_2p0K_11$koel_dis, ann_dis_2p0K_12$koel_dis,
                               ann_dis_2p0K_13$koel_dis)

ann_mel_2p0K_base <- apply(ann_mel_2p0K_base_all, 1, mea_na)
ann_mel_2p0K_coch <- apply(ann_mel_2p0K_coch_all, 1, mea_na)
ann_mel_2p0K_koel <- apply(ann_mel_2p0K_koel_all, 1, mea_na)
ann_lpr_2p0K_base <- apply(ann_lpr_2p0K_base_all, 1, mea_na)
ann_lpr_2p0K_coch <- apply(ann_lpr_2p0K_coch_all, 1, mea_na)
ann_lpr_2p0K_koel <- apply(ann_lpr_2p0K_koel_all, 1, mea_na)
ann_fra_2p0K_base <- apply(ann_fra_2p0K_base_all, 1, mea_na)
ann_fra_2p0K_coch <- apply(ann_fra_2p0K_coch_all, 1, mea_na)
ann_fra_2p0K_koel <- apply(ann_fra_2p0K_koel_all, 1, mea_na)
ann_ele_2p0K_base <- apply(ann_ele_2p0K_base_all, 1, mea_na)
ann_ele_2p0K_coch <- apply(ann_ele_2p0K_coch_all, 1, mea_na)
ann_ele_2p0K_koel <- apply(ann_ele_2p0K_koel_all, 1, mea_na)

ann_prt_2p0K_base <- apply(ann_prt_2p0K_base_all, 1, mea_na)
ann_prt_2p0K_coch <- apply(ann_prt_2p0K_coch_all, 1, mea_na)
ann_prt_2p0K_koel <- apply(ann_prt_2p0K_koel_all, 1, mea_na)
ann_prl_2p0K_base <- apply(ann_prl_2p0K_base_all, 1, mea_na)
ann_prl_2p0K_coch <- apply(ann_prl_2p0K_coch_all, 1, mea_na)
ann_prl_2p0K_koel <- apply(ann_prl_2p0K_koel_all, 1, mea_na)
ann_eff_2p0K_base <- apply(ann_eff_2p0K_base_all, 1, mea_na)
ann_eff_2p0K_coch <- apply(ann_eff_2p0K_coch_all, 1, mea_na)
ann_eff_2p0K_koel <- apply(ann_eff_2p0K_koel_all, 1, mea_na)

ann_dis_2p0K_base <- apply(ann_dis_2p0K_base_all, 1, mea_na)
ann_dis_2p0K_coch <- apply(ann_dis_2p0K_coch_all, 1, mea_na)
ann_dis_2p0K_koel <- apply(ann_dis_2p0K_koel_all, 1, mea_na)

#3.0K warming level
ann_flu_3p0K_1 <-  annu_flu(flux_3p0K_1)
ann_flu_3p0K_2 <-  annu_flu(flux_3p0K_2)
ann_flu_3p0K_3 <-  annu_flu(flux_3p0K_3)
ann_flu_3p0K_4 <-  annu_flu(flux_3p0K_4)
ann_flu_3p0K_5 <-  annu_flu(flux_3p0K_5)
ann_flu_3p0K_6 <-  annu_flu(flux_3p0K_6)
ann_flu_3p0K_7 <-  annu_flu(flux_3p0K_7)
ann_flu_3p0K_8 <-  annu_flu(flux_3p0K_8)

ann_pre_3p0K_1 <-  annu_pre(prec_3p0K_1)
ann_pre_3p0K_2 <-  annu_pre(prec_3p0K_2)
ann_pre_3p0K_3 <-  annu_pre(prec_3p0K_3)
ann_pre_3p0K_4 <-  annu_pre(prec_3p0K_4)
ann_pre_3p0K_5 <-  annu_pre(prec_3p0K_5)
ann_pre_3p0K_6 <-  annu_pre(prec_3p0K_6)
ann_pre_3p0K_7 <-  annu_pre(prec_3p0K_7)
ann_pre_3p0K_8 <-  annu_pre(prec_3p0K_8)

ann_dis_3p0K_1 <-  annu_dis(disc_3p0K_1)
ann_dis_3p0K_2 <-  annu_dis(disc_3p0K_2)
ann_dis_3p0K_3 <-  annu_dis(disc_3p0K_3)
ann_dis_3p0K_4 <-  annu_dis(disc_3p0K_4)
ann_dis_3p0K_5 <-  annu_dis(disc_3p0K_5)
ann_dis_3p0K_6 <-  annu_dis(disc_3p0K_6)
ann_dis_3p0K_7 <-  annu_dis(disc_3p0K_7)
ann_dis_3p0K_8 <-  annu_dis(disc_3p0K_8)

ann_mel_3p0K_base_all <- cbind(ann_flu_3p0K_1$mel_base_mea,  ann_flu_3p0K_2$mel_base_mea,  ann_flu_3p0K_3$mel_base_mea,
                               ann_flu_3p0K_4$mel_base_mea,  ann_flu_3p0K_5$mel_base_mea,  ann_flu_3p0K_6$mel_base_mea,
                               ann_flu_3p0K_7$mel_base_mea,  ann_flu_3p0K_8$mel_base_mea)
ann_mel_3p0K_coch_all <- cbind(ann_flu_3p0K_1$mel_coch_mea,  ann_flu_3p0K_2$mel_coch_mea,  ann_flu_3p0K_3$mel_coch_mea,
                               ann_flu_3p0K_4$mel_coch_mea,  ann_flu_3p0K_5$mel_coch_mea,  ann_flu_3p0K_6$mel_coch_mea,
                               ann_flu_3p0K_7$mel_coch_mea,  ann_flu_3p0K_8$mel_coch_mea)
ann_mel_3p0K_koel_all <- cbind(ann_flu_3p0K_1$mel_koel_mea,  ann_flu_3p0K_2$mel_koel_mea,  ann_flu_3p0K_3$mel_koel_mea,
                               ann_flu_3p0K_4$mel_koel_mea,  ann_flu_3p0K_5$mel_koel_mea,  ann_flu_3p0K_6$mel_koel_mea,
                               ann_flu_3p0K_7$mel_koel_mea,  ann_flu_3p0K_8$mel_koel_mea)
ann_lpr_3p0K_base_all <- cbind(ann_flu_3p0K_1$lpr_base_mea,  ann_flu_3p0K_2$lpr_base_mea,  ann_flu_3p0K_3$lpr_base_mea,
                               ann_flu_3p0K_4$lpr_base_mea,  ann_flu_3p0K_5$lpr_base_mea,  ann_flu_3p0K_6$lpr_base_mea,
                               ann_flu_3p0K_7$lpr_base_mea,  ann_flu_3p0K_8$lpr_base_mea)
ann_lpr_3p0K_coch_all <- cbind(ann_flu_3p0K_1$lpr_coch_mea,  ann_flu_3p0K_2$lpr_coch_mea,  ann_flu_3p0K_3$lpr_coch_mea,
                               ann_flu_3p0K_4$lpr_coch_mea,  ann_flu_3p0K_5$lpr_coch_mea,  ann_flu_3p0K_6$lpr_coch_mea,
                               ann_flu_3p0K_7$lpr_coch_mea,  ann_flu_3p0K_8$lpr_coch_mea)
ann_lpr_3p0K_koel_all <- cbind(ann_flu_3p0K_1$lpr_koel_mea,  ann_flu_3p0K_2$lpr_koel_mea,  ann_flu_3p0K_3$lpr_koel_mea,
                               ann_flu_3p0K_4$lpr_koel_mea,  ann_flu_3p0K_5$lpr_koel_mea,  ann_flu_3p0K_6$lpr_koel_mea,
                               ann_flu_3p0K_7$lpr_koel_mea,  ann_flu_3p0K_8$lpr_koel_mea)
ann_fra_3p0K_base_all <- cbind(ann_flu_3p0K_1$fra_base_mea,  ann_flu_3p0K_2$fra_base_mea,  ann_flu_3p0K_3$fra_base_mea,
                               ann_flu_3p0K_4$fra_base_mea,  ann_flu_3p0K_5$fra_base_mea,  ann_flu_3p0K_6$fra_base_mea,
                               ann_flu_3p0K_7$fra_base_mea,  ann_flu_3p0K_8$fra_base_mea)
ann_fra_3p0K_coch_all <- cbind(ann_flu_3p0K_1$fra_coch_mea,  ann_flu_3p0K_2$fra_coch_mea,  ann_flu_3p0K_3$fra_coch_mea,
                               ann_flu_3p0K_4$fra_coch_mea,  ann_flu_3p0K_5$fra_coch_mea,  ann_flu_3p0K_6$fra_coch_mea,
                               ann_flu_3p0K_7$fra_coch_mea,  ann_flu_3p0K_8$fra_coch_mea)
ann_fra_3p0K_koel_all <- cbind(ann_flu_3p0K_1$fra_koel_mea,  ann_flu_3p0K_2$fra_koel_mea,  ann_flu_3p0K_3$fra_koel_mea,
                               ann_flu_3p0K_4$fra_koel_mea,  ann_flu_3p0K_5$fra_koel_mea,  ann_flu_3p0K_6$fra_koel_mea,
                               ann_flu_3p0K_7$fra_koel_mea,  ann_flu_3p0K_8$fra_koel_mea)
ann_ele_3p0K_base_all <- cbind(ann_flu_3p0K_1$ele_base_mea,  ann_flu_3p0K_2$ele_base_mea,  ann_flu_3p0K_3$ele_base_mea,
                               ann_flu_3p0K_4$ele_base_mea,  ann_flu_3p0K_5$ele_base_mea,  ann_flu_3p0K_6$ele_base_mea,
                               ann_flu_3p0K_7$ele_base_mea,  ann_flu_3p0K_8$ele_base_mea)
ann_ele_3p0K_coch_all <- cbind(ann_flu_3p0K_1$ele_coch_mea,  ann_flu_3p0K_2$ele_coch_mea,  ann_flu_3p0K_3$ele_coch_mea,
                               ann_flu_3p0K_4$ele_coch_mea,  ann_flu_3p0K_5$ele_coch_mea,  ann_flu_3p0K_6$ele_coch_mea,
                               ann_flu_3p0K_7$ele_coch_mea,  ann_flu_3p0K_8$ele_coch_mea)
ann_ele_3p0K_koel_all <- cbind(ann_flu_3p0K_1$ele_koel_mea,  ann_flu_3p0K_2$ele_koel_mea,  ann_flu_3p0K_3$ele_koel_mea,
                               ann_flu_3p0K_4$ele_koel_mea,  ann_flu_3p0K_5$ele_koel_mea,  ann_flu_3p0K_6$ele_koel_mea,
                               ann_flu_3p0K_7$ele_koel_mea,  ann_flu_3p0K_8$ele_koel_mea)

ann_prt_3p0K_base_all <- cbind(ann_pre_3p0K_1$base_pre_tot,  ann_pre_3p0K_2$base_pre_tot,  ann_pre_3p0K_3$base_pre_tot,
                               ann_pre_3p0K_4$base_pre_tot,  ann_pre_3p0K_5$base_pre_tot,  ann_pre_3p0K_6$base_pre_tot,
                               ann_pre_3p0K_7$base_pre_tot,  ann_pre_3p0K_8$base_pre_tot)
ann_prt_3p0K_coch_all <- cbind(ann_pre_3p0K_1$coch_pre_tot,  ann_pre_3p0K_2$coch_pre_tot,  ann_pre_3p0K_3$coch_pre_tot,
                               ann_pre_3p0K_4$coch_pre_tot,  ann_pre_3p0K_5$coch_pre_tot,  ann_pre_3p0K_6$coch_pre_tot,
                               ann_pre_3p0K_7$coch_pre_tot,  ann_pre_3p0K_8$coch_pre_tot)
ann_prt_3p0K_koel_all <- cbind(ann_pre_3p0K_1$koel_pre_tot,  ann_pre_3p0K_2$koel_pre_tot,  ann_pre_3p0K_3$koel_pre_tot,
                               ann_pre_3p0K_4$koel_pre_tot,  ann_pre_3p0K_5$koel_pre_tot,  ann_pre_3p0K_6$koel_pre_tot,
                               ann_pre_3p0K_7$koel_pre_tot,  ann_pre_3p0K_8$koel_pre_tot)
ann_prl_3p0K_base_all <- cbind(ann_pre_3p0K_1$base_pre_liq,  ann_pre_3p0K_2$base_pre_liq,  ann_pre_3p0K_3$base_pre_liq,
                               ann_pre_3p0K_4$base_pre_liq,  ann_pre_3p0K_5$base_pre_liq,  ann_pre_3p0K_6$base_pre_liq,
                               ann_pre_3p0K_7$base_pre_liq,  ann_pre_3p0K_8$base_pre_liq)
ann_prl_3p0K_coch_all <- cbind(ann_pre_3p0K_1$coch_pre_liq,  ann_pre_3p0K_2$coch_pre_liq,  ann_pre_3p0K_3$coch_pre_liq,
                               ann_pre_3p0K_4$coch_pre_liq,  ann_pre_3p0K_5$coch_pre_liq,  ann_pre_3p0K_6$coch_pre_liq,
                               ann_pre_3p0K_7$coch_pre_liq,  ann_pre_3p0K_8$coch_pre_liq)
ann_prl_3p0K_koel_all <- cbind(ann_pre_3p0K_1$koel_pre_liq,  ann_pre_3p0K_2$koel_pre_liq,  ann_pre_3p0K_3$koel_pre_liq,
                               ann_pre_3p0K_4$koel_pre_liq,  ann_pre_3p0K_5$koel_pre_liq,  ann_pre_3p0K_6$koel_pre_liq,
                               ann_pre_3p0K_7$koel_pre_liq,  ann_pre_3p0K_8$koel_pre_liq)
ann_eff_3p0K_base_all <- cbind(ann_pre_3p0K_1$base_pro_eff,  ann_pre_3p0K_2$base_pro_eff,  ann_pre_3p0K_3$base_pro_eff,
                               ann_pre_3p0K_4$base_pro_eff,  ann_pre_3p0K_5$base_pro_eff,  ann_pre_3p0K_6$base_pro_eff,
                               ann_pre_3p0K_7$base_pro_eff,  ann_pre_3p0K_8$base_pro_eff)
ann_eff_3p0K_coch_all <- cbind(ann_pre_3p0K_1$coch_pro_eff,  ann_pre_3p0K_2$coch_pro_eff,  ann_pre_3p0K_3$coch_pro_eff,
                               ann_pre_3p0K_4$coch_pro_eff,  ann_pre_3p0K_5$coch_pro_eff,  ann_pre_3p0K_6$coch_pro_eff,
                               ann_pre_3p0K_7$coch_pro_eff,  ann_pre_3p0K_8$coch_pro_eff)
ann_eff_3p0K_koel_all <- cbind(ann_pre_3p0K_1$koel_pro_eff,  ann_pre_3p0K_2$koel_pro_eff,  ann_pre_3p0K_3$koel_pro_eff,
                               ann_pre_3p0K_4$koel_pro_eff,  ann_pre_3p0K_5$koel_pro_eff,  ann_pre_3p0K_6$koel_pro_eff,
                               ann_pre_3p0K_7$koel_pro_eff,  ann_pre_3p0K_8$koel_pro_eff)

ann_dis_3p0K_base_all <- cbind(ann_dis_3p0K_1$base_dis,  ann_dis_3p0K_2$base_dis,  ann_dis_3p0K_3$base_dis,
                               ann_dis_3p0K_4$base_dis,  ann_dis_3p0K_5$base_dis,  ann_dis_3p0K_6$base_dis,
                               ann_dis_3p0K_7$base_dis,  ann_dis_3p0K_8$base_dis)
ann_dis_3p0K_coch_all <- cbind(ann_dis_3p0K_1$coch_dis,  ann_dis_3p0K_2$coch_dis,  ann_dis_3p0K_3$coch_dis,
                               ann_dis_3p0K_4$coch_dis,  ann_dis_3p0K_5$coch_dis,  ann_dis_3p0K_6$coch_dis,
                               ann_dis_3p0K_7$coch_dis,  ann_dis_3p0K_8$coch_dis)
ann_dis_3p0K_koel_all <- cbind(ann_dis_3p0K_1$koel_dis,  ann_dis_3p0K_2$koel_dis,  ann_dis_3p0K_3$koel_dis,
                               ann_dis_3p0K_4$koel_dis,  ann_dis_3p0K_5$koel_dis,  ann_dis_3p0K_6$koel_dis,
                               ann_dis_3p0K_7$koel_dis,  ann_dis_3p0K_8$koel_dis)

ann_mel_3p0K_base <- apply(ann_mel_3p0K_base_all, 1, mea_na)
ann_mel_3p0K_coch <- apply(ann_mel_3p0K_coch_all, 1, mea_na)
ann_mel_3p0K_koel <- apply(ann_mel_3p0K_koel_all, 1, mea_na)
ann_lpr_3p0K_base <- apply(ann_lpr_3p0K_base_all, 1, mea_na)
ann_lpr_3p0K_coch <- apply(ann_lpr_3p0K_coch_all, 1, mea_na)
ann_lpr_3p0K_koel <- apply(ann_lpr_3p0K_koel_all, 1, mea_na)
ann_fra_3p0K_base <- apply(ann_fra_3p0K_base_all, 1, mea_na)
ann_fra_3p0K_coch <- apply(ann_fra_3p0K_coch_all, 1, mea_na)
ann_fra_3p0K_koel <- apply(ann_fra_3p0K_koel_all, 1, mea_na)
ann_ele_3p0K_base <- apply(ann_ele_3p0K_base_all, 1, mea_na)
ann_ele_3p0K_coch <- apply(ann_ele_3p0K_coch_all, 1, mea_na)
ann_ele_3p0K_koel <- apply(ann_ele_3p0K_koel_all, 1, mea_na)

ann_prt_3p0K_base <- apply(ann_prt_3p0K_base_all, 1, mea_na)
ann_prt_3p0K_coch <- apply(ann_prt_3p0K_coch_all, 1, mea_na)
ann_prt_3p0K_koel <- apply(ann_prt_3p0K_koel_all, 1, mea_na)
ann_prl_3p0K_base <- apply(ann_prl_3p0K_base_all, 1, mea_na)
ann_prl_3p0K_coch <- apply(ann_prl_3p0K_coch_all, 1, mea_na)
ann_prl_3p0K_koel <- apply(ann_prl_3p0K_koel_all, 1, mea_na)
ann_eff_3p0K_base <- apply(ann_eff_3p0K_base_all, 1, mea_na)
ann_eff_3p0K_coch <- apply(ann_eff_3p0K_coch_all, 1, mea_na)
ann_eff_3p0K_koel <- apply(ann_eff_3p0K_koel_all, 1, mea_na)

ann_dis_3p0K_base <- apply(ann_dis_3p0K_base_all, 1, mea_na)
ann_dis_3p0K_coch <- apply(ann_dis_3p0K_coch_all, 1, mea_na)
ann_dis_3p0K_koel <- apply(ann_dis_3p0K_koel_all, 1, mea_na)


#Functiont to plot annual cycles
ann_cycl_2 <- function(data_hist_all, data_1p5K_all, data_2p0K_all, data_3p0K_all,
                       data_hist_mea, data_1p5K_mea, data_2p0K_mea, data_3p0K_mea, 
                       main = "", do_legend = F,
                       y_lab = "", do_y_lab = F){
  
  col_hist <- "steelblue4"
  col_1p5K <- "grey25"
  col_2p0K <- "orange3"
  col_3p0K <- "darkred"
  alpha_points <- 0.06
  alpha_points_mea <- 0.6
  cex_points <- 0.6
  cex_x_label <- 1.6
  cex_main <- 1.7
  
  ylims <- c(min_na(c(data_hist_all, data_1p5K_all, data_2p0K_all, data_3p0K_all)),
             max_na(c(data_hist_all, data_1p5K_all, data_2p0K_all, data_3p0K_all)))
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  plot(data_hist_all[, 1], type = "n", col = "blue3", ylim = ylims, ylab = "", xlab = "", axes = F)
  abline(v = x_axis_tic, lwd = 0.8, col = "grey55", lty = "dashed")
  for (i in 1:ncol(data_hist_all)) {
    points(data_hist_all[, i], col = scales::alpha(col_hist, alpha_points), pch = 19, cex = cex_points)  
  }
  for (i in 1:ncol(data_1p5K_all)) {
    points(data_1p5K_all[, i], col = scales::alpha(col_1p5K, alpha_points), pch = 19, cex = cex_points)
  }
  for (i in 1:ncol(data_2p0K_all)) {
    points(data_2p0K_all[, i], col = scales::alpha(col_2p0K, alpha_points), pch = 19, cex = cex_points)  
  }
  for (i in 1:ncol(data_3p0K_all)) {
    points(data_3p0K_all[, i], col = scales::alpha(col_3p0K, alpha_points), pch = 19, cex = cex_points)  
  }
  points(data_hist_mea, col = scales::alpha(col_hist, alpha_points_mea), pch = 19, cex = cex_points)
  points(data_1p5K_mea, col = scales::alpha(col_1p5K, alpha_points_mea), pch = 19, cex = cex_points)
  points(data_2p0K_mea, col = scales::alpha(col_2p0K, alpha_points_mea), pch = 19, cex = cex_points)
  points(data_3p0K_mea, col = scales::alpha(col_3p0K, alpha_points_mea), pch = 19, cex = cex_points)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.05)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.40, 0), cex.axis = cex_x_label)#plot labels
  axis(2, mgp=c(3, 0.29, 0), tck = -0.015, cex.axis = 1.8)
  mtext(main, side = 3, line = 0.3, cex = cex_main, adj = 0.0)
  if(do_y_lab){
    mtext(y_lab, side = 3, line = 1.8, cex = 1.3)
  }
  
  if(do_legend){
    legend("topright", c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
           col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.6,
           box.lwd = 0.2, box.col = "black", bg = "white")
  }
  box()
  
}




#Plot annual cycles export
pdf(paste0(bas_dir,"res_figs/ann_cyc_fut3.pdf"), width = 16, height = 6)
# tiff(paste0(bas_dir, "res_figs/ann_cyc_fut3",".tiff"), width = 10.0, height = 6.0,
#      units = "in", res = 300)

par(family = "serif")
par(mar = c(1.8, 3.0, 1.0, 1.0))

layout(matrix(c(7, 7, 7, 7,
                8, 1, 2, 3,
                8, 4, 5, 6),
              3, 4, byrow = T), widths=c(0.1, 1, 1, 1), heights=c(0.15, rep(1, 2)))
# layout.show(n = 8)


#Protective effect
ann_cycl_2(ann_eff_hist_base_all, ann_eff_1p5K_base_all, 
           ann_eff_2p0K_base_all, ann_eff_3p0K_base_all,
           ann_eff_hist_base, ann_eff_1p5K_base, ann_eff_2p0K_base, ann_eff_3p0K_base,
           y_lab = "", do_y_lab = T, do_legend = T) 

#Melt fraction
ann_cycl_2(ann_fra_hist_base_all, ann_fra_1p5K_base_all, 
           ann_fra_2p0K_base_all, ann_fra_3p0K_base_all,
           ann_fra_hist_base, ann_fra_1p5K_base, ann_fra_2p0K_base, ann_fra_3p0K_base,
           y_lab = "", do_y_lab = T) 

#Snowmelt elevation
ann_cycl_2(ann_ele_hist_base_all, ann_ele_1p5K_base_all, 
           ann_ele_2p0K_base_all, ann_ele_3p0K_base_all,
           ann_ele_hist_base, ann_ele_1p5K_base, ann_ele_2p0K_base, ann_ele_3p0K_base,
           y_lab = "", do_y_lab = T) 

ann_cycl_2(ann_eff_hist_coch_all, ann_eff_1p5K_coch_all, 
           ann_eff_2p0K_coch_all, ann_eff_3p0K_coch_all,
           ann_eff_hist_coch, ann_eff_1p5K_coch, ann_eff_2p0K_coch, ann_eff_3p0K_coch)

ann_cycl_2(ann_fra_hist_coch_all, ann_fra_1p5K_coch_all, 
           ann_fra_2p0K_coch_all, ann_fra_3p0K_coch_all,
           ann_fra_hist_coch, ann_fra_1p5K_coch, ann_fra_2p0K_coch, ann_fra_3p0K_coch)

ann_cycl_2(ann_ele_hist_coch_all, ann_ele_1p5K_coch_all, 
           ann_ele_2p0K_coch_all, ann_ele_3p0K_coch_all,
           ann_ele_hist_coch, ann_ele_1p5K_coch, ann_ele_2p0K_coch, ann_ele_3p0K_coch)

cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Protective effect [-]",  side = 3, line = -2.4, cex = cex_header, adj = 0.14, outer = T)
mtext("2. Melt fraction [-]",  side = 3, line = -2.4, cex = cex_header, adj = 0.52, outer = T)
mtext("3. Melt elevation [m]",  side = 3, line = -2.4, cex = cex_header, adj = 0.90, outer = T)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 2, line = -2.9, cex = cex_header+0.2, adj = 0.83)
mtext("b) Cochem",
      side = 2, line = -2.9, cex = cex_header+0.2, adj = 0.17)


dev.off()








pdf(paste0(bas_dir,"res_figs/ann_cyc_fut2.pdf"), width = 16, height = 16)

par(family = "serif")
par(mar = c(1.5, 3.0, 1.0, 1.0))

layout(matrix(c(rep(22, 4),
                23, 1, 2, 3,
                23, 4, 5, 6,
                23, 7, 8, 9,
                23, 10, 11, 12,
                23, 13, 14, 15,
                23, 16, 17, 18,
                23, 19, 20, 21),
              8, 4, byrow = T), widths=c(0.15, 1, 1, 1), heights=c(0.2, rep(1, 7)))
# layout.show(n = 23)

#Discharge
ann_cycl_2(ann_dis_hist_base_all, ann_dis_1p5K_base_all, 
           ann_dis_2p0K_base_all, ann_dis_3p0K_base_all,
           ann_dis_hist_base, ann_dis_1p5K_base, ann_dis_2p0K_base, ann_dis_3p0K_base,
           y_lab = expression(paste("Discharge [m"^"3", "s"^"-1","]")), 
           do_legend = T, do_y_lab = T)

ann_cycl_2(ann_dis_hist_coch_all, ann_dis_1p5K_coch_all, 
           ann_dis_2p0K_coch_all, ann_dis_3p0K_coch_all,
           ann_dis_hist_coch, ann_dis_1p5K_coch, ann_dis_2p0K_coch, ann_dis_3p0K_coch,
           y_lab = expression(paste("Discharge [m"^"3", "s"^"-1","]")),
           do_legend = T)

ann_cycl_2(ann_dis_hist_koel_all, ann_dis_1p5K_koel_all, 
           ann_dis_2p0K_koel_all, ann_dis_3p0K_koel_all,
           ann_dis_hist_koel, ann_dis_1p5K_koel, ann_dis_2p0K_koel, ann_dis_3p0K_koel,
           y_lab = expression(paste("Discharge [m"^"3", "s"^"-1","]")),
           do_legend = T)

#Precipitation total
ann_cycl_2(ann_prt_hist_base_all, ann_prt_1p5K_base_all, 
           ann_prt_2p0K_base_all, ann_prt_3p0K_base_all,
           ann_prt_hist_base, ann_prt_1p5K_base, ann_prt_2p0K_base, ann_prt_3p0K_base,
           y_lab = expression(paste("Precip. [mm]")), 
           do_y_lab = T)

ann_cycl_2(ann_prt_hist_coch_all, ann_prt_1p5K_coch_all, 
           ann_prt_2p0K_coch_all, ann_prt_3p0K_coch_all,
           ann_prt_hist_coch, ann_prt_1p5K_coch, ann_prt_2p0K_coch, ann_prt_3p0K_coch)

ann_cycl_2(ann_prt_hist_koel_all, ann_prt_1p5K_koel_all, 
           ann_prt_2p0K_koel_all, ann_prt_3p0K_koel_all,
           ann_prt_hist_koel, ann_prt_1p5K_koel, ann_prt_2p0K_koel, ann_prt_3p0K_koel)

#Precipitation total
ann_cycl_2(ann_prl_hist_base_all, ann_prl_1p5K_base_all, 
           ann_prl_2p0K_base_all, ann_prl_3p0K_base_all,
           ann_prl_hist_base, ann_prl_1p5K_base, ann_prl_2p0K_base, ann_prl_3p0K_base,
           y_lab = expression(paste("Precip. [mm]")), 
           do_y_lab = T)

ann_cycl_2(ann_prl_hist_coch_all, ann_prl_1p5K_coch_all, 
           ann_prl_2p0K_coch_all, ann_prl_3p0K_coch_all,
           ann_prl_hist_coch, ann_prl_1p5K_coch, ann_prl_2p0K_coch, ann_prl_3p0K_coch)

ann_cycl_2(ann_prl_hist_koel_all, ann_prl_1p5K_koel_all, 
           ann_prl_2p0K_koel_all, ann_prl_3p0K_koel_all,
           ann_prl_hist_koel, ann_prl_1p5K_koel, ann_prl_2p0K_koel, ann_prl_3p0K_koel)

#Protective effect
ann_cycl_2(ann_eff_hist_base_all, ann_eff_1p5K_base_all, 
           ann_eff_2p0K_base_all, ann_eff_3p0K_base_all,
           ann_eff_hist_base, ann_eff_1p5K_base, ann_eff_2p0K_base, ann_eff_3p0K_base,
           y_lab = expression(paste("Prec. solid/total [-]")), do_y_lab = T) 

ann_cycl_2(ann_eff_hist_coch_all, ann_eff_1p5K_coch_all, 
           ann_eff_2p0K_coch_all, ann_eff_3p0K_coch_all,
           ann_eff_hist_coch, ann_eff_1p5K_coch, ann_eff_2p0K_coch, ann_eff_3p0K_coch)

ann_cycl_2(ann_eff_hist_koel_all, ann_eff_1p5K_koel_all, 
           ann_eff_2p0K_koel_all, ann_eff_3p0K_koel_all,
           ann_eff_hist_koel, ann_eff_1p5K_koel, ann_eff_2p0K_koel, ann_eff_3p0K_koel)

#Snowmelt
ann_cycl_2(ann_mel_hist_base_all, ann_mel_1p5K_base_all, 
           ann_mel_2p0K_base_all, ann_mel_3p0K_base_all,
           ann_mel_hist_base, ann_mel_1p5K_base, ann_mel_2p0K_base, ann_mel_3p0K_base,
           y_lab = expression(paste("Snowmelt [mm]")), do_y_lab = T) 

ann_cycl_2(ann_mel_hist_coch_all, ann_mel_1p5K_coch_all, 
           ann_mel_2p0K_coch_all, ann_mel_3p0K_coch_all,
           ann_mel_hist_coch, ann_mel_1p5K_coch, ann_mel_2p0K_coch, ann_mel_3p0K_coch)

ann_cycl_2(ann_mel_hist_koel_all, ann_mel_1p5K_koel_all, 
           ann_mel_2p0K_koel_all, ann_mel_3p0K_koel_all,
           ann_mel_hist_koel, ann_mel_1p5K_koel, ann_mel_2p0K_koel, ann_mel_3p0K_koel)

#Runoff fraction
ann_cycl_2(ann_fra_hist_base_all, ann_fra_1p5K_base_all, 
           ann_fra_2p0K_base_all, ann_fra_3p0K_base_all,
           ann_fra_hist_base, ann_fra_1p5K_base, ann_fra_2p0K_base, ann_fra_3p0K_base,
           y_lab = expression(paste("Snowmelt fraction [-]")), do_y_lab = T) 

ann_cycl_2(ann_fra_hist_coch_all, ann_fra_1p5K_coch_all, 
           ann_fra_2p0K_coch_all, ann_fra_3p0K_coch_all,
           ann_fra_hist_coch, ann_fra_1p5K_coch, ann_fra_2p0K_coch, ann_fra_3p0K_coch)

ann_cycl_2(ann_fra_hist_koel_all, ann_fra_1p5K_koel_all, 
           ann_fra_2p0K_koel_all, ann_fra_3p0K_koel_all,
           ann_fra_hist_koel, ann_fra_1p5K_koel, ann_fra_2p0K_koel, ann_fra_3p0K_koel)

#Snowmelt elevation
ann_cycl_2(ann_ele_hist_base_all, ann_ele_1p5K_base_all, 
           ann_ele_2p0K_base_all, ann_ele_3p0K_base_all,
           ann_ele_hist_base, ann_ele_1p5K_base, ann_ele_2p0K_base, ann_ele_3p0K_base,
           y_lab = expression(paste("Elevation [-]")), do_y_lab = T) 

ann_cycl_2(ann_ele_hist_coch_all, ann_ele_1p5K_coch_all, 
           ann_ele_2p0K_coch_all, ann_ele_3p0K_coch_all,
           ann_ele_hist_coch, ann_ele_1p5K_coch, ann_ele_2p0K_coch, ann_ele_3p0K_coch)

ann_cycl_2(ann_ele_hist_koel_all, ann_ele_1p5K_koel_all, 
           ann_ele_2p0K_koel_all, ann_ele_3p0K_koel_all,
           ann_ele_hist_koel, ann_ele_1p5K_koel, ann_ele_2p0K_koel, ann_ele_3p0K_koel)

cex_header <- 1.7
par(mar = c(0,0,0,0))

#Gauging station

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.191)
mtext("b) Cochem",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.525)
mtext("c) Cologne",
      side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.875)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Disc.",  side = 2, line = -2.2, cex = cex_header, adj = 0.933, outer = T)
mtext("2. Prec. total",  side = 2, line = -2.2, cex = cex_header, adj = 0.805, outer = T)
mtext("3. Prec. liquid",  side = 2, line = -2.2, cex = cex_header, adj = 0.643, outer = T)
mtext("4. Protect. effect",  side = 2, line = -2.2, cex = cex_header, adj = 0.483, outer = T)
mtext("5. Snowmelt",  side = 2, line = -2.2, cex = cex_header, adj = 0.325, outer = T)
mtext("6. Melt fract.",  side = 2, line = -2.2, cex = cex_header, adj = 0.180, outer = T)
mtext("7. Melt elevat.",  side = 2, line = -2.2, cex = cex_header, adj = 0.024, outer = T)

dev.off()




#snow_future----

#Snow cover duration simulations

scd_from_nc <- function(gcm_model, delta_t, rcp){
  
  #select nc_file
  nc_path_sel <- nc_file_paths[which(grepl(gcm_model, nc_file_paths) &
                                       grepl(rcp, nc_file_paths))]
  nc_file_sel <- nc_open(paste0(nc_path_sel, "output/mHM_Fluxes_States.nc"))
  
  #path to precipitation input
  # nc_path_sel_prec <- gsub("output", "input/meteo", nc_path_sel)
  # pr_file_sel <- nc_open(paste0(nc_path_sel, "input/pre.nc"))
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  
  #if simulation time frame does not entirely cover warming period
  if(date_sel[1] > date[1]){
    sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
    count_date <- length(date_sel)
  }else{
    sta_date_ind <- 1
    count_date <- length(date_sel) - which(format(date_sel) == date[1]) + 1
  }
  
  #snowpack
  snow_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
  
  sd2sc <- function(val_in, sc_thr = 3){
    
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
  
  for(i in 1:count_date){
    
    # print(i)
    
    scd_sim_sing <- sapply(c(snow_cube[, , i]), sd2sc)
    
    if(i == 1){
      scd_sim <-  scd_sim_sing
    }else{
      scd_sim <-  scd_sim + scd_sim_sing
    }
    
  }
  
  #annual average snow cover duration
  scd_sim_ann <- scd_sim / round(count_date / 365)
  
  return(scd_sim_ann)
  
}

#historical
tic()
scd_hist_1 <- scd_from_nc("GFDL-ESM2M",     "historical", "historical")
scd_hist_2 <- scd_from_nc("HadGEM2-ES",     "historical", "historical")
scd_hist_3 <- scd_from_nc("IPSL-CM5A-LR",   "historical", "historical")
scd_hist_4 <- scd_from_nc("MIROC-ESM-CHEM", "historical", "historical")
scd_hist_5 <- scd_from_nc("NorESM1-M",      "historical", "historical")
toc()

#1.5K warming level
tic()
scd_1p5K_1  <- scd_from_nc("HadGEM2-ES",     "1p5", "2p6")
scd_1p5K_2  <- scd_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
scd_1p5K_3  <- scd_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
scd_1p5K_4  <- scd_from_nc("NorESM1-M",      "1p5", "2p6")
scd_1p5K_5  <- scd_from_nc("GFDL-ESM2M",     "1p5", "6p0")
scd_1p5K_6  <- scd_from_nc("HadGEM2-ES",     "1p5", "6p0")
scd_1p5K_7  <- scd_from_nc("IPSL-CM5A-LR",   "1p5", "6p0")
scd_1p5K_8  <- scd_from_nc("MIROC-ESM-CHEM", "1p5", "6p0")
scd_1p5K_9  <- scd_from_nc("NorESM1-M",      "1p5", "6p0")
scd_1p5K_10 <- scd_from_nc("GFDL-ESM2M",     "1p5", "8p5")
scd_1p5K_11 <- scd_from_nc("HadGEM2-ES",     "1p5", "8p5")
scd_1p5K_12 <- scd_from_nc("IPSL-CM5A-LR",   "1p5", "8p5")
scd_1p5K_13 <- scd_from_nc("MIROC-ESM-CHEM", "1p5", "8p5")
scd_1p5K_14 <- scd_from_nc("NorESM1-M",      "1p5", "8p5")
toc()

#2K warming level
tic()
scd_2p0K_1  <- scd_from_nc("HadGEM2-ES",     "2p0", "2p6")
scd_2p0K_2  <- scd_from_nc("IPSL-CM5A-LR",   "2p0", "2p6")
scd_2p0K_3  <- scd_from_nc("MIROC-ESM-CHEM", "2p0", "2p6")
scd_2p0K_4  <- scd_from_nc("GFDL-ESM2M",     "2p0", "6p0")
scd_2p0K_5  <- scd_from_nc("HadGEM2-ES",     "2p0", "6p0")
scd_2p0K_6  <- scd_from_nc("IPSL-CM5A-LR",   "2p0", "6p0")
scd_2p0K_7  <- scd_from_nc("MIROC-ESM-CHEM", "2p0", "6p0")
scd_2p0K_8  <- scd_from_nc("NorESM1-M",      "2p0", "6p0")
scd_2p0K_9  <- scd_from_nc("GFDL-ESM2M",     "2p0", "8p5")
scd_2p0K_10 <- scd_from_nc("HadGEM2-ES",     "2p0", "8p5")
scd_2p0K_11 <- scd_from_nc("IPSL-CM5A-LR",   "2p0", "8p5")
scd_2p0K_12 <- scd_from_nc("MIROC-ESM-CHEM", "2p0", "8p5")
scd_2p0K_13 <- scd_from_nc("NorESM1-M",      "2p0", "8p5")
toc()

#3K warming level
tic()
scd_3p0K_1 <- scd_from_nc("HadGEM2-ES",     "3p0", "6p0")
scd_3p0K_2 <- scd_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
scd_3p0K_3 <- scd_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
scd_3p0K_4 <- scd_from_nc("GFDL-ESM2M",     "3p0", "8p5")
scd_3p0K_5 <- scd_from_nc("HadGEM2-ES",     "3p0", "8p5")
scd_3p0K_6 <- scd_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
scd_3p0K_7 <- scd_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
scd_3p0K_8 <- scd_from_nc("NorESM1-M",      "3p0", "8p5")
toc()

scd_hist_all <- cbind(scd_hist_1,  scd_hist_2,  scd_hist_3,  scd_hist_4, scd_hist_5)
scd_1p5K_all <- cbind(scd_1p5K_1,  scd_1p5K_2,  scd_1p5K_3,  scd_1p5K_4, scd_1p5K_5,
                      scd_1p5K_6,  scd_1p5K_7,  scd_1p5K_8,  scd_1p5K_9, scd_1p5K_10, 
                      scd_1p5K_11, scd_1p5K_12, scd_1p5K_13, scd_1p5K_14)
scd_2p0K_all <- cbind(scd_2p0K_1,  scd_2p0K_2,  scd_2p0K_3,  scd_2p0K_4, scd_2p0K_5,
                      scd_2p0K_6,  scd_2p0K_7,  scd_2p0K_8,  scd_2p0K_9, scd_2p0K_10, 
                      scd_2p0K_11, scd_2p0K_12, scd_2p0K_13)
scd_3p0K_all <- cbind(scd_3p0K_1,  scd_3p0K_2,  scd_3p0K_3,  scd_3p0K_4, scd_3p0K_5,
                      scd_3p0K_6,  scd_3p0K_7,  scd_3p0K_8)

scd_hist <- apply(scd_hist_all, 1, mea_na)
scd_1p5K <- apply(scd_1p5K_all, 1, mea_na)
scd_2p0K <- apply(scd_2p0K_all, 1, mea_na)
scd_3p0K <- apply(scd_3p0K_all, 1, mea_na)

#get cells in basin
nc_dummy <- nc_open(paste0(run_dir, "output/GFDL-ESM2M/historical/output/mRM_Fluxes_States.nc"))
lon <- ncdf4::ncvar_get(nc_dummy, varid = "lon")
lat <- ncdf4::ncvar_get(nc_dummy, varid = "lat")

#spatial grid points from lat/lon
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)

n_cores <- 5 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

cols_spat_sim_hist <- foreach(i = 1:length(scd_hist), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_hist[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}
cols_spat_sim_1p5K <- foreach(i = 1:length(scd_1p5K), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_1p5K[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}
cols_spat_sim_2p0K <- foreach(i = 1:length(scd_2p0K), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_2p0K[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}
cols_spat_sim_3p0K <- foreach(i = 1:length(scd_3p0K), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_3p0K[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}

scd_dif_1 <- (scd_2p0K - scd_hist)
scd_dif_2 <- (scd_3p0K - scd_hist)
scd_dif_3 <- (scd_3p0K - scd_2p0K)
scd_lims <- range(c(scd_dif_1, scd_dif_2, scd_dif_3), na.rm = T)

cols_spat_dif_1 <- foreach(i = 1:length(scd_dif_1), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif_1[i], 
          dat_ref = scd_lims[1]:scd_lims[2],
          do_log = F,
          do_bicol = T)
  
}
cols_spat_dif_2 <- foreach(i = 1:length(scd_dif_2), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif_2[i], 
          dat_ref = scd_lims,
          do_log = F,
          do_bicol = T)
  
}
cols_spat_dif_3 <- foreach(i = 1:length(scd_dif_3), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif_3[i], 
          dat_ref = scd_lims,
          do_log = F,
          do_bicol = T)
  
}

stopCluster(my_clust)


pdf(paste0(bas_dir, "res_figs/scd_maps_sim",".pdf"), width = 16, height = 2*4.2)

layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4,  rep(5, 7),  6,
                rep(7, 7), 8, rep(9, 7), 10, rep(11, 7), 12),
              2, 24, byrow = T), widths=c(), heights=c())
# layout.show(n = 12)

par(family = "serif")
cex_pch <- 0.60

#Historic
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_sim_hist, cex = cex_pch)
mtext("a) Historic", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, 365, length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_hist), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#2.0K
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_sim_2p0K, cex = cex_pch)
mtext("b) 2.0K", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, 365, length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_2p0K), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#3.0K
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_sim_3p0K, cex = cex_pch)
mtext("c) 3.0K", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, 365, length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_3p0K), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Difference Historic to 2.0K
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_dif_1, cex = cex_pch)
mtext("d) Diff. 2.0K - Hist.", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(scd_lims[1], -scd_lims[1], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif_1), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Difference Historic to 3.0K
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_dif_2, cex = cex_pch)
mtext("e) Diff. 3.0K - Hist.", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(scd_lims[1], -scd_lims[1], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif_2), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Difference 2.0K to 3.0K
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[, 1], grid_points_cube_84@coords[, 2], pch = 15,
       col = cols_spat_dif_3, cex = cex_pch)
mtext("f) Diff. 3.0K - 2.0K", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(scd_lims[1], -scd_lims[1], length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif_3), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

dev.off()



#snow_vali----

#Make cluster for parallel computing
n_cores <- 5 #number of cores used for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

#Snow cover duration MODIS

scf_dlr_dir <- "D:/nrc_user/rottler/SCF_data/snow_dlr/SnowPack_DLR.tar/SnowPack_DLR/" 

file_names <- dir(path = scf_dlr_dir, recursive = T)
# unique(nchar(file_names))
# file_names[which(nchar(file_names) %in% c(36, 32, 10))]
file_names <- file_names[which(nchar(file_names) == nchar(file_names[1]))]
scf_file <- raster(paste0(scf_dlr_dir , file_names[1]))

basin_lobi_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/lobith_catch.shp")
basin_lobi_raw_84 <- spTransform(basin_lobi_buf, CRS = crswgs84)

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
#aggregate to mHM resolution
scf_buf_aggr <- aggregate(scf_buf, fact = 10, fun = mean, na.rm = TRUE)
plot(scf_buf, col = viridis(200, direction = -1))
plot(scf_buf_aggr, col = viridis(200, direction = -1))

#Snow cover duration EOBS simulations
nc_flux_file <- paste0(run_dir, "output/EOBS/output/mHM_Fluxes_States.nc")
nc_flux <- nc_open(nc_flux_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

count_date <- length(date)

#Fluxes and states
snow_cube <- ncvar_get(nc_flux, start = c(1, 1, 1), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")

sd2sc <- function(val_in, sc_thr = 3){
  
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

date_scd <- scf_date #MODIS data coverage

date_scd_ind <- which(date %in% date_scd)

cells_sel <- which(!is.na(c(snow_cube[, , date_scd_ind[1]])))

for(i in 1:length(date_scd_ind)){
  
  print(i)
  
  scd_sim_sing <- sapply(c(snow_cube[, , date_scd_ind[i]][cells_sel]), sd2sc)
  
  if(i == 1){
    scd_eob_sum <-  scd_sim_sing
  }else{
    scd_eob_sum <-  scd_eob_sum + scd_sim_sing
  }
  
}

scd_eob <- scd_eob_sum / round(length(scf_date)/365)

#Get values MODIS grid points simulated
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)
scd_dlr_sum <- raster::extract(scf_buf_aggr, grid_points_cube_84[cells_sel])
scd_dlr <- scd_dlr_sum / (round(length(scf_date)/365))

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
cols_spat_sim <- foreach(i = 1:length(scd_eob), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_eob[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}

#Values to colors difference
scd_dif <- (scd_eob - scd_dlr)
cols_spat_dif <- foreach(i = 1:length(scd_dif), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif[i], 
          dat_ref = scd_dif,
          do_log = F,
          do_bicol = T)
  
}

#Values to colors observations
cols_spat_obs <- foreach(i = 1:length(scd_dlr), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dlr[i],
          dat_ref = 0:365,
          do_bicol = F)
  
}

stopCluster(my_clust)

pdf(paste0(bas_dir, "res_figs/scd_maps_vali",".pdf"), width = 16, height = 4.2)

#Plot maps
layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4, rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())
# layout.show(n = 7)

par(family = "serif")
cex_pch <- 0.60

#Map Simulations
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_lobi_raw_84, border = alpha("black", alpha = 0.0))
points(grid_points_cube_84@coords[cells_sel, 1], grid_points_cube_84@coords[cells_sel, 2], 
       pch = 15, col = cols_spat_sim, cex = cex_pch)
# plot(basin_base, add =T, lwd = 1.5)
mtext("a) EOBS simulations", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, 365, length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_eob), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
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
my_bre <- seq(0, 365, length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dlr_ann), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

dev.off()











#snow_tower----

#Snow cover duration EOBS simulations
nc_flux_file <- paste0(run_dir, "output/EOBS/output/mHM_Fluxes_States.nc")
nc_flux <- nc_open(nc_flux_file)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_flux, varid = "lon")
lat <- ncdf4::ncvar_get(nc_flux, varid = "lat")
date <- as.Date(as.character(nc.get.time.series(nc_flux, time.dim.name = "time")))

count_date <- length(date)

#Fluxes and states
snow_cube <- ncvar_get(nc_flux, start = c(1, 1, 1), 
                       count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")

snow_max <- apply(snow_cube, c(1, 2), max_na) #[mm]

snow_max_c <- c(snow_max)

sn_tow_1 <- which(snow_max_c > 2500)

#Plot snow towers
# sn_tow_1
st_sel_ind <- c(11324, 11411, 11412)

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

stow_ts_1 <- ncvar_get(nc_flux, start = c(sel_x_1, sel_y_1, 1), 
                       count = c(1, 1, count_date), varid = "snowpack")

stow_ts_2 <- ncvar_get(nc_flux, start = c(sel_x_2, sel_y_2, 1), 
                       count = c(1, 1, count_date), varid = "snowpack")

stow_ts_3 <- ncvar_get(nc_flux, start = c(sel_x_3, sel_y_3, 1), 
                       count = c(1, 1, count_date), varid = "snowpack")


pdf(paste0(bas_dir,"res_figs/snow_tower.pdf"), width = 12, height = 6)

par(family = "serif")
par(mfrow = c(3, 1))
par(mar = c(2, 4, 3, 1))

plot(date, stow_ts_1, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[1]], 3), "  Lon: ", round(c(lon)[st_sel_ind[1]], 3)),
      side = 3, line = 0.2, cex = 1.2)

plot(date, stow_ts_2, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[2]], 3), "  Lon: ", round(c(lon)[st_sel_ind[2]], 3)),
      side = 3, line = 0.2, cex = 1.2)

plot(date, stow_ts_3, type = "l", ylab = "", cex.axis = 1.2, col = "darkblue", lwd = 1.5)
mtext("SWE [mm]", side = 2, line = 2.8)
mtext(paste0("Lat: ", round(c(lat)[st_sel_ind[3]], 3), "  Lon: ", round(c(lon)[st_sel_ind[3]], 3)),
      side = 3, line = 0.2, cex = 1.2)

dev.off()


#number of snow towers in GCM-driven simulations
number_towers <- function(gcm_model, delta_t, rcp, ncores = 5){
  
  #select nc_file
  nc_path_sel <- nc_file_paths[which(grepl(gcm_model, nc_file_paths) &
                                       grepl(rcp, nc_file_paths))]
  nc_file_sel <- nc_open(paste0(nc_path_sel, "output/mHM_Fluxes_States.nc"))
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  
  #if simulation time frame does not entirely cover warming period
  if(date_sel[1] > date[1]){
    sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
    count_date <- length(date_sel)
  }else{
    sta_date_ind <- 1
    count_date <- length(date_sel) - which(format(date_sel) == date[1]) + 1
  }
  
  #snowpack
  snow_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
  
  snow_max <- apply(snow_cube, c(1, 2), max_na) #[mm]
  
  snow_max_c <- c(snow_max)
  
  sn_tow <- length(which(snow_max_c > 2500))
  
  return(sn_tow)
  
}

#historical
tic()
stow_hist_1 <- number_towers("GFDL-ESM2M",     "historical", "historical")
stow_hist_2 <- number_towers("HadGEM2-ES",     "historical", "historical")
stow_hist_3 <- number_towers("IPSL-CM5A-LR",   "historical", "historical")
stow_hist_4 <- number_towers("MIROC-ESM-CHEM", "historical", "historical")
stow_hist_5 <- number_towers("NorESM1-M",      "historical", "historical")
toc()


#1.5K warming level
tic()
stow_1p5K_1  <- number_towers("HadGEM2-ES",     "1p5", "2p6")
stow_1p5K_2  <- number_towers("IPSL-CM5A-LR",   "1p5", "2p6")
stow_1p5K_3  <- number_towers("MIROC-ESM-CHEM", "1p5", "2p6")
stow_1p5K_4  <- number_towers("NorESM1-M",      "1p5", "2p6")
stow_1p5K_5  <- number_towers("GFDL-ESM2M",     "1p5", "6p0")
stow_1p5K_6  <- number_towers("HadGEM2-ES",     "1p5", "6p0")
stow_1p5K_7  <- number_towers("IPSL-CM5A-LR",   "1p5", "6p0")
stow_1p5K_8  <- number_towers("MIROC-ESM-CHEM", "1p5", "6p0")
stow_1p5K_9  <- number_towers("NorESM1-M",      "1p5", "6p0")
stow_1p5K_10 <- number_towers("GFDL-ESM2M",     "1p5", "8p5")
stow_1p5K_11 <- number_towers("HadGEM2-ES",     "1p5", "8p5")
stow_1p5K_12 <- number_towers("IPSL-CM5A-LR",   "1p5", "8p5")
stow_1p5K_13 <- number_towers("MIROC-ESM-CHEM", "1p5", "8p5")
stow_1p5K_14 <- number_towers("NorESM1-M",      "1p5", "8p5")
toc()


#2K warming level
tic()
stow_2p0K_1  <- number_towers("HadGEM2-ES",     "2p0", "2p6")
stow_2p0K_2  <- number_towers("IPSL-CM5A-LR",   "2p0", "2p6")
stow_2p0K_3  <- number_towers("MIROC-ESM-CHEM", "2p0", "2p6")
stow_2p0K_4  <- number_towers("GFDL-ESM2M",     "2p0", "6p0")
stow_2p0K_5  <- number_towers("HadGEM2-ES",     "2p0", "6p0")
stow_2p0K_6  <- number_towers("IPSL-CM5A-LR",   "2p0", "6p0")
stow_2p0K_7  <- number_towers("MIROC-ESM-CHEM", "2p0", "6p0")
stow_2p0K_8  <- number_towers("NorESM1-M",      "2p0", "6p0")
stow_2p0K_9  <- number_towers("GFDL-ESM2M",     "2p0", "8p5")
stow_2p0K_10 <- number_towers("HadGEM2-ES",     "2p0", "8p5")
stow_2p0K_11 <- number_towers("IPSL-CM5A-LR",   "2p0", "8p5")
stow_2p0K_12 <- number_towers("MIROC-ESM-CHEM", "2p0", "8p5")
stow_2p0K_13 <- number_towers("NorESM1-M",      "2p0", "8p5")
toc()


#3K warming level
tic()
stow_3p0K_1 <- number_towers("HadGEM2-ES",     "3p0", "6p0")
stow_3p0K_2 <- number_towers("IPSL-CM5A-LR",   "3p0", "6p0")
stow_3p0K_3 <- number_towers("MIROC-ESM-CHEM", "3p0", "6p0")
stow_3p0K_4 <- number_towers("GFDL-ESM2M",     "3p0", "8p5")
stow_3p0K_5 <- number_towers("HadGEM2-ES",     "3p0", "8p5")
stow_3p0K_6 <- number_towers("IPSL-CM5A-LR",   "3p0", "8p5")
stow_3p0K_7 <- number_towers("MIROC-ESM-CHEM", "3p0", "8p5")
stow_3p0K_8 <- number_towers("NorESM1-M",      "3p0", "8p5")
toc()

#month_max----

f_max_month <- function(data_in){
  
  date_in <- data_in$date #get date from input table
  data_in$date <- NULL #remove date column
  sta_year <- as.numeric(format(date_in[1], '%Y'))
  end_year <- as.numeric(format(date_in[length(date_in)], '%Y'))
  
  max_out <- NULL
  for(s in 1:ncol(data_in)){
    
    data_day <- meltimr::ord_day(data_in = data_in[, s], 
                                 date = date_in, 
                                 start_y = sta_year,
                                 end_y = end_year)
    
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
    month_cols <- list(oct_cols, nov_cols, dec_cols, jan_cols, feb_cols, mar_cols, apr_cols, 
                       may_cols, jun_cols, jul_cols, aug_cols, sep_cols)
    
    max_sta <- NULL
    for(y in 1:nrow(data_day)){
      
      max_mon <- NULL
      for(i in 1:length(month_cols)){
        
        max_mon <- c(max_mon, max_na(data_day[y, month_cols[[i]]]))
        
      }  
      
      max_sta <- rbind(max_sta, max_mon)
      
    }  
    
    max_out <- cbind(max_out, max_sta)
    
  }
  
  return(max_out)
  
}

#Montly maxima: Historical
dmax_mon_hist_1 <- f_max_month(disc_hist_1) ; fmax_mon_hist_1 <- f_max_month(flux_hist_1)
dmax_mon_hist_2 <- f_max_month(disc_hist_2) ; fmax_mon_hist_2 <- f_max_month(flux_hist_2)
dmax_mon_hist_3 <- f_max_month(disc_hist_3) ; fmax_mon_hist_3 <- f_max_month(flux_hist_3)
dmax_mon_hist_4 <- f_max_month(disc_hist_4) ; fmax_mon_hist_4 <- f_max_month(flux_hist_4)
dmax_mon_hist_5 <- f_max_month(disc_hist_5) ; fmax_mon_hist_5 <- f_max_month(flux_hist_5)

pmax_mon_hist_1 <- f_max_month(prec_hist_1)
pmax_mon_hist_2 <- f_max_month(prec_hist_2)
pmax_mon_hist_3 <- f_max_month(prec_hist_3)
pmax_mon_hist_4 <- f_max_month(prec_hist_4)
pmax_mon_hist_5 <- f_max_month(prec_hist_5)

#Montly maxima: 1.5K warming
dmax_mon_1p5K_1 <- f_max_month(disc_1p5K_1) ; fmax_mon_1p5K_1 <- f_max_month(flux_1p5K_1) 
dmax_mon_1p5K_2 <- f_max_month(disc_1p5K_2) ; fmax_mon_1p5K_2 <- f_max_month(flux_1p5K_2)
dmax_mon_1p5K_3 <- f_max_month(disc_1p5K_3) ; fmax_mon_1p5K_3 <- f_max_month(flux_1p5K_3)
dmax_mon_1p5K_4 <- f_max_month(disc_1p5K_4) ; fmax_mon_1p5K_4 <- f_max_month(flux_1p5K_4)
dmax_mon_1p5K_5 <- f_max_month(disc_1p5K_5) ; fmax_mon_1p5K_5 <- f_max_month(flux_1p5K_5)
dmax_mon_1p5K_6 <- f_max_month(disc_1p5K_6) ; fmax_mon_1p5K_6 <- f_max_month(flux_1p5K_6)
dmax_mon_1p5K_7 <- f_max_month(disc_1p5K_7) ; fmax_mon_1p5K_7 <- f_max_month(flux_1p5K_7)
dmax_mon_1p5K_8 <- f_max_month(disc_1p5K_8) ; fmax_mon_1p5K_8 <- f_max_month(flux_1p5K_8)
dmax_mon_1p5K_9 <- f_max_month(disc_1p5K_9) ; fmax_mon_1p5K_9 <- f_max_month(flux_1p5K_9)
dmax_mon_1p5K_10 <- f_max_month(disc_1p5K_10) ; fmax_mon_1p5K_10 <- f_max_month(flux_1p5K_10)
dmax_mon_1p5K_11 <- f_max_month(disc_1p5K_11) ; fmax_mon_1p5K_11 <- f_max_month(flux_1p5K_11)
dmax_mon_1p5K_12 <- f_max_month(disc_1p5K_12) ; fmax_mon_1p5K_12 <- f_max_month(flux_1p5K_12)
dmax_mon_1p5K_13 <- f_max_month(disc_1p5K_13) ; fmax_mon_1p5K_13 <- f_max_month(flux_1p5K_13)
dmax_mon_1p5K_14 <- f_max_month(disc_1p5K_14) ; fmax_mon_1p5K_14 <- f_max_month(flux_1p5K_14)

pmax_mon_1p5K_1 <- f_max_month(prec_1p5K_1)
pmax_mon_1p5K_2 <- f_max_month(prec_1p5K_2)
pmax_mon_1p5K_3 <- f_max_month(prec_1p5K_3)
pmax_mon_1p5K_4 <- f_max_month(prec_1p5K_4)
pmax_mon_1p5K_5 <- f_max_month(prec_1p5K_5)
pmax_mon_1p5K_6 <- f_max_month(prec_1p5K_6)
pmax_mon_1p5K_7 <- f_max_month(prec_1p5K_7)
pmax_mon_1p5K_8 <- f_max_month(prec_1p5K_8)
pmax_mon_1p5K_9 <- f_max_month(prec_1p5K_9) 
pmax_mon_1p5K_10 <- f_max_month(prec_1p5K_10)
pmax_mon_1p5K_11 <- f_max_month(prec_1p5K_11)
pmax_mon_1p5K_12 <- f_max_month(prec_1p5K_12)
pmax_mon_1p5K_13 <- f_max_month(prec_1p5K_13)
pmax_mon_1p5K_14 <- f_max_month(prec_1p5K_14)

#Monthly maxima: 2.0K warming
dmax_mon_2p0K_1 <- f_max_month(disc_2p0K_1) ; fmax_mon_2p0K_1 <- f_max_month(flux_2p0K_1)
dmax_mon_2p0K_2 <- f_max_month(disc_2p0K_2) ; fmax_mon_2p0K_2 <- f_max_month(flux_2p0K_2)
dmax_mon_2p0K_3 <- f_max_month(disc_2p0K_3) ; fmax_mon_2p0K_3 <- f_max_month(flux_2p0K_3)
dmax_mon_2p0K_4 <- f_max_month(disc_2p0K_4) ; fmax_mon_2p0K_4 <- f_max_month(flux_2p0K_4)
dmax_mon_2p0K_5 <- f_max_month(disc_2p0K_5) ; fmax_mon_2p0K_5 <- f_max_month(flux_2p0K_5)
dmax_mon_2p0K_6 <- f_max_month(disc_2p0K_6) ; fmax_mon_2p0K_6 <- f_max_month(flux_2p0K_6)
dmax_mon_2p0K_7 <- f_max_month(disc_2p0K_7) ; fmax_mon_2p0K_7 <- f_max_month(flux_2p0K_7)
dmax_mon_2p0K_8 <- f_max_month(disc_2p0K_8) ; fmax_mon_2p0K_8 <- f_max_month(flux_2p0K_8)
dmax_mon_2p0K_9 <- f_max_month(disc_2p0K_9) ; fmax_mon_2p0K_9 <- f_max_month(flux_2p0K_9)
dmax_mon_2p0K_10 <- f_max_month(disc_2p0K_10) ; fmax_mon_2p0K_10 <- f_max_month(flux_2p0K_10)
dmax_mon_2p0K_11 <- f_max_month(disc_2p0K_11) ; fmax_mon_2p0K_11 <- f_max_month(flux_2p0K_11)
dmax_mon_2p0K_12 <- f_max_month(disc_2p0K_12) ; fmax_mon_2p0K_12 <- f_max_month(flux_2p0K_12)
dmax_mon_2p0K_13 <- f_max_month(disc_2p0K_13) ; fmax_mon_2p0K_13 <- f_max_month(flux_2p0K_13)

pmax_mon_2p0K_1 <- f_max_month(prec_2p0K_1)
pmax_mon_2p0K_2 <- f_max_month(prec_2p0K_2) 
pmax_mon_2p0K_3 <- f_max_month(prec_2p0K_3) 
pmax_mon_2p0K_4 <- f_max_month(prec_2p0K_4)
pmax_mon_2p0K_5 <- f_max_month(prec_2p0K_5) 
pmax_mon_2p0K_6 <- f_max_month(prec_2p0K_6) 
pmax_mon_2p0K_7 <- f_max_month(prec_2p0K_7) 
pmax_mon_2p0K_8 <- f_max_month(prec_2p0K_8)
pmax_mon_2p0K_9 <- f_max_month(prec_2p0K_9)
pmax_mon_2p0K_10 <- f_max_month(prec_2p0K_10) 
pmax_mon_2p0K_11 <- f_max_month(prec_2p0K_11) 
pmax_mon_2p0K_12 <- f_max_month(prec_2p0K_12) 
pmax_mon_2p0K_13 <- f_max_month(prec_2p0K_13)

#Monthly maxima: 3.0K warming
dmax_mon_3p0K_1 <- f_max_month(disc_3p0K_1) ; fmax_mon_3p0K_1 <- f_max_month(flux_3p0K_1)
dmax_mon_3p0K_2 <- f_max_month(disc_3p0K_2) ; fmax_mon_3p0K_2 <- f_max_month(flux_3p0K_2)
dmax_mon_3p0K_3 <- f_max_month(disc_3p0K_3) ; fmax_mon_3p0K_3 <- f_max_month(flux_3p0K_3)
dmax_mon_3p0K_4 <- f_max_month(disc_3p0K_4) ; fmax_mon_3p0K_4 <- f_max_month(flux_3p0K_4)
dmax_mon_3p0K_5 <- f_max_month(disc_3p0K_5) ; fmax_mon_3p0K_5 <- f_max_month(flux_3p0K_5)
dmax_mon_3p0K_6 <- f_max_month(disc_3p0K_6) ; fmax_mon_3p0K_6 <- f_max_month(flux_3p0K_6)
dmax_mon_3p0K_7 <- f_max_month(disc_3p0K_7) ; fmax_mon_3p0K_7 <- f_max_month(flux_3p0K_7)
dmax_mon_3p0K_8 <- f_max_month(disc_3p0K_8) ; fmax_mon_3p0K_8 <- f_max_month(flux_3p0K_8)

pmax_mon_3p0K_1 <- f_max_month(prec_3p0K_1)
pmax_mon_3p0K_2 <- f_max_month(prec_3p0K_2)
pmax_mon_3p0K_3 <- f_max_month(prec_3p0K_3)
pmax_mon_3p0K_4 <- f_max_month(prec_3p0K_4)
pmax_mon_3p0K_5 <- f_max_month(prec_3p0K_5)
pmax_mon_3p0K_6 <- f_max_month(prec_3p0K_6)
pmax_mon_3p0K_7 <- f_max_month(prec_3p0K_7)
pmax_mon_3p0K_8 <- f_max_month(prec_3p0K_8)

#Put together historical
max_dis_mon_hist_koel <- rbind(dmax_mon_hist_1[, 1:12], dmax_mon_hist_2[, 1:12], dmax_mon_hist_3[, 1:12],
                               dmax_mon_hist_4[, 1:12], dmax_mon_hist_5[, 1:12])
max_dis_mon_hist_coch <- rbind(dmax_mon_hist_1[, 13:24], dmax_mon_hist_2[, 13:24], dmax_mon_hist_3[, 13:24],
                               dmax_mon_hist_4[, 13:24], dmax_mon_hist_5[, 13:24])
max_dis_mon_hist_base <- rbind(dmax_mon_hist_1[, 25:36], dmax_mon_hist_2[, 25:36], dmax_mon_hist_3[, 25:36],
                               dmax_mon_hist_4[, 25:36], dmax_mon_hist_5[, 25:36])

max_mel_mon_hist_koel <- rbind(fmax_mon_hist_1[, 25:36], fmax_mon_hist_2[, 25:36], fmax_mon_hist_3[, 25:36],
                               fmax_mon_hist_4[, 25:36], fmax_mon_hist_5[, 25:36])
max_mel_mon_hist_coch <- rbind(fmax_mon_hist_1[, 13:24], fmax_mon_hist_2[, 13:24], fmax_mon_hist_3[, 13:24],
                               fmax_mon_hist_4[, 13:24], fmax_mon_hist_5[, 13:24])
max_mel_mon_hist_base <- rbind(fmax_mon_hist_1[, 1:12], fmax_mon_hist_2[, 1:12], fmax_mon_hist_3[, 1:12],
                               fmax_mon_hist_4[, 1:12], fmax_mon_hist_5[, 1:12])

max_prt_mon_hist_koel <- rbind(pmax_mon_hist_1[, 25:36], pmax_mon_hist_2[, 25:36], pmax_mon_hist_3[, 25:36],
                               pmax_mon_hist_4[, 25:36], pmax_mon_hist_5[, 25:36])
max_prt_mon_hist_coch <- rbind(pmax_mon_hist_1[, 13:24], pmax_mon_hist_2[, 13:24], pmax_mon_hist_3[, 13:24],
                               pmax_mon_hist_4[, 13:24], pmax_mon_hist_5[, 13:24])
max_prt_mon_hist_base <- rbind(pmax_mon_hist_1[, 1:12], pmax_mon_hist_2[, 1:12], pmax_mon_hist_3[, 1:12],
                               pmax_mon_hist_4[, 1:12], pmax_mon_hist_5[, 1:12])

max_prl_mon_hist_koel <- rbind(pmax_mon_hist_1[, 61:72], pmax_mon_hist_2[, 61:72], pmax_mon_hist_3[, 61:72],
                               pmax_mon_hist_4[, 61:72], pmax_mon_hist_5[, 61:72])
max_prl_mon_hist_coch <- rbind(pmax_mon_hist_1[, 49:60], pmax_mon_hist_2[, 49:60], pmax_mon_hist_3[, 49:60],
                               pmax_mon_hist_4[, 49:60], pmax_mon_hist_5[, 49:60])
max_prl_mon_hist_base <- rbind(pmax_mon_hist_1[, 37:48], pmax_mon_hist_2[, 37:48], pmax_mon_hist_3[, 37:48],
                               pmax_mon_hist_4[, 37:48], pmax_mon_hist_5[, 37:48])

#Put together 1.5K warming level
max_dis_mon_1p5K_koel <- rbind(dmax_mon_1p5K_1[, 1:12], dmax_mon_1p5K_2[, 1:12], dmax_mon_1p5K_3[, 1:12],
                               dmax_mon_1p5K_4[, 1:12], dmax_mon_1p5K_5[, 1:12], dmax_mon_1p5K_6[, 1:12],
                               dmax_mon_1p5K_7[, 1:12], dmax_mon_1p5K_8[, 1:12], dmax_mon_1p5K_9[, 1:12],
                               dmax_mon_1p5K_10[, 1:12], dmax_mon_1p5K_11[, 1:12], dmax_mon_1p5K_12[, 1:12],
                               dmax_mon_1p5K_13[, 1:12], dmax_mon_1p5K_14[, 1:12])
max_dis_mon_1p5K_coch <- rbind(dmax_mon_1p5K_1[, 13:24], dmax_mon_1p5K_2[, 13:24], dmax_mon_1p5K_3[, 13:24],
                               dmax_mon_1p5K_4[, 13:24], dmax_mon_1p5K_5[, 13:24], dmax_mon_1p5K_6[, 13:24],
                               dmax_mon_1p5K_7[, 13:24], dmax_mon_1p5K_8[, 13:24], dmax_mon_1p5K_9[, 13:24],
                               dmax_mon_1p5K_10[, 13:24], dmax_mon_1p5K_11[, 13:24], dmax_mon_1p5K_12[, 13:24],
                               dmax_mon_1p5K_13[, 13:24], dmax_mon_1p5K_14[, 13:24])
max_dis_mon_1p5K_base <- rbind(dmax_mon_1p5K_1[, 25:36], dmax_mon_1p5K_2[, 25:36], dmax_mon_1p5K_3[, 25:36],
                               dmax_mon_1p5K_4[, 25:36], dmax_mon_1p5K_5[, 25:36], dmax_mon_1p5K_6[, 25:36],
                               dmax_mon_1p5K_7[, 25:36], dmax_mon_1p5K_8[, 25:36], dmax_mon_1p5K_9[, 25:36],
                               dmax_mon_1p5K_10[, 25:36], dmax_mon_1p5K_11[, 25:36], dmax_mon_1p5K_12[, 25:36],
                               dmax_mon_1p5K_13[, 25:36], dmax_mon_1p5K_14[, 25:36])

max_mel_mon_1p5K_koel <- rbind(fmax_mon_1p5K_1[, 25:36], fmax_mon_1p5K_2[, 25:36], fmax_mon_1p5K_3[, 25:36],
                               fmax_mon_1p5K_4[, 25:36], fmax_mon_1p5K_5[, 25:36], fmax_mon_1p5K_6[, 25:36],
                               fmax_mon_1p5K_7[, 25:36], fmax_mon_1p5K_8[, 25:36], fmax_mon_1p5K_9[, 25:36],
                               fmax_mon_1p5K_10[, 25:36], fmax_mon_1p5K_11[, 25:36], fmax_mon_1p5K_12[, 25:36],
                               fmax_mon_1p5K_13[, 25:36], fmax_mon_1p5K_14[, 25:36])
max_mel_mon_1p5K_coch <- rbind(fmax_mon_1p5K_1[, 13:24], fmax_mon_1p5K_2[, 13:24], fmax_mon_1p5K_3[, 13:24],
                               fmax_mon_1p5K_4[, 13:24], fmax_mon_1p5K_5[, 13:24], fmax_mon_1p5K_6[, 13:24],
                               fmax_mon_1p5K_7[, 13:24], fmax_mon_1p5K_8[, 13:24], fmax_mon_1p5K_9[, 13:24],
                               fmax_mon_1p5K_10[, 13:24], fmax_mon_1p5K_11[, 13:24], fmax_mon_1p5K_12[, 13:24],
                               fmax_mon_1p5K_13[, 13:24], fmax_mon_1p5K_14[, 13:24])
max_mel_mon_1p5K_base <- rbind(fmax_mon_1p5K_1[, 1:12], fmax_mon_1p5K_2[, 1:12], fmax_mon_1p5K_3[, 1:12],
                               fmax_mon_1p5K_4[, 1:12], fmax_mon_1p5K_5[, 1:12], fmax_mon_1p5K_6[, 1:12],
                               fmax_mon_1p5K_7[, 1:12], fmax_mon_1p5K_8[, 1:12], fmax_mon_1p5K_9[, 1:12],
                               fmax_mon_1p5K_10[, 1:12], fmax_mon_1p5K_11[, 1:12], fmax_mon_1p5K_12[, 1:12],
                               fmax_mon_1p5K_13[, 1:12], fmax_mon_1p5K_14[, 1:12])

max_prt_mon_1p5K_koel <- rbind(pmax_mon_1p5K_1[, 25:36], pmax_mon_1p5K_2[, 25:36], pmax_mon_1p5K_3[, 25:36],
                               pmax_mon_1p5K_4[, 25:36], pmax_mon_1p5K_5[, 25:36], pmax_mon_1p5K_6[, 25:36],
                               pmax_mon_1p5K_7[, 25:36], pmax_mon_1p5K_8[, 25:36], pmax_mon_1p5K_9[, 25:36],
                               pmax_mon_1p5K_10[, 25:36], pmax_mon_1p5K_11[, 25:36], pmax_mon_1p5K_12[, 25:36],
                               pmax_mon_1p5K_13[, 25:36], pmax_mon_1p5K_14[, 25:36])
max_prt_mon_1p5K_coch <- rbind(pmax_mon_1p5K_1[, 13:24], pmax_mon_1p5K_2[, 13:24], pmax_mon_1p5K_3[, 13:24],
                               pmax_mon_1p5K_4[, 13:24], pmax_mon_1p5K_5[, 13:24], pmax_mon_1p5K_6[, 13:24],
                               pmax_mon_1p5K_7[, 13:24], pmax_mon_1p5K_8[, 13:24], pmax_mon_1p5K_9[, 13:24],
                               pmax_mon_1p5K_10[, 13:24], pmax_mon_1p5K_11[, 13:24], pmax_mon_1p5K_12[, 13:24],
                               pmax_mon_1p5K_13[, 13:24], pmax_mon_1p5K_14[, 13:24])
max_prt_mon_1p5K_base <- rbind(pmax_mon_1p5K_1[, 1:12], pmax_mon_1p5K_2[, 1:12], pmax_mon_1p5K_3[, 1:12],
                               pmax_mon_1p5K_4[, 1:12], pmax_mon_1p5K_5[, 1:12], pmax_mon_1p5K_6[, 1:12],
                               pmax_mon_1p5K_7[, 1:12], pmax_mon_1p5K_8[, 1:12], pmax_mon_1p5K_9[, 1:12],
                               pmax_mon_1p5K_10[, 1:12], pmax_mon_1p5K_11[, 1:12], pmax_mon_1p5K_12[, 1:12],
                               pmax_mon_1p5K_13[, 1:12], pmax_mon_1p5K_14[, 1:12])

max_prl_mon_1p5K_koel <- rbind(pmax_mon_1p5K_1[, 61:72], pmax_mon_1p5K_2[, 61:72], pmax_mon_1p5K_3[, 61:72],
                               pmax_mon_1p5K_4[, 61:72], pmax_mon_1p5K_5[, 61:72], pmax_mon_1p5K_6[, 61:72],
                               pmax_mon_1p5K_7[, 61:72], pmax_mon_1p5K_8[, 61:72], pmax_mon_1p5K_9[, 61:72],
                               pmax_mon_1p5K_10[, 61:72], pmax_mon_1p5K_11[, 61:72], pmax_mon_1p5K_12[, 61:72],
                               pmax_mon_1p5K_13[, 61:72], pmax_mon_1p5K_14[, 61:72])
max_prl_mon_1p5K_coch <- rbind(pmax_mon_1p5K_1[, 49:60], pmax_mon_1p5K_2[, 49:60], pmax_mon_1p5K_3[, 49:60],
                               pmax_mon_1p5K_4[, 49:60], pmax_mon_1p5K_5[, 49:60], pmax_mon_1p5K_6[, 49:60],
                               pmax_mon_1p5K_7[, 49:60], pmax_mon_1p5K_8[, 49:60], pmax_mon_1p5K_9[, 49:60],
                               pmax_mon_1p5K_10[, 49:60], pmax_mon_1p5K_11[, 49:60], pmax_mon_1p5K_12[, 49:60],
                               pmax_mon_1p5K_13[, 49:60], pmax_mon_1p5K_14[, 49:60])
max_prl_mon_1p5K_base <- rbind(pmax_mon_1p5K_1[, 37:48], pmax_mon_1p5K_2[, 37:48], pmax_mon_1p5K_3[, 37:48],
                               pmax_mon_1p5K_4[, 37:48], pmax_mon_1p5K_5[, 37:48], pmax_mon_1p5K_6[, 37:48],
                               pmax_mon_1p5K_7[, 37:48], pmax_mon_1p5K_8[, 37:48], pmax_mon_1p5K_9[, 37:48],
                               pmax_mon_1p5K_10[, 37:48], pmax_mon_1p5K_11[, 37:48], pmax_mon_1p5K_12[, 37:48],
                               pmax_mon_1p5K_13[, 37:48], pmax_mon_1p5K_14[, 37:48])

#Put together 2.0K warming level
max_dis_mon_2p0K_koel <- rbind(dmax_mon_2p0K_1[, 1:12], dmax_mon_2p0K_2[, 1:12], dmax_mon_2p0K_3[, 1:12],
                               dmax_mon_2p0K_4[, 1:12], dmax_mon_2p0K_5[, 1:12], dmax_mon_2p0K_6[, 1:12],
                               dmax_mon_2p0K_7[, 1:12], dmax_mon_2p0K_8[, 1:12], dmax_mon_2p0K_9[, 1:12],
                               dmax_mon_2p0K_10[, 1:12], dmax_mon_2p0K_11[, 1:12], dmax_mon_2p0K_12[, 1:12],
                               dmax_mon_2p0K_13[, 1:12])
max_dis_mon_2p0K_coch <- rbind(dmax_mon_2p0K_1[, 13:24], dmax_mon_2p0K_2[, 13:24], dmax_mon_2p0K_3[, 13:24],
                               dmax_mon_2p0K_4[, 13:24], dmax_mon_2p0K_5[, 13:24], dmax_mon_2p0K_6[, 13:24],
                               dmax_mon_2p0K_7[, 13:24], dmax_mon_2p0K_8[, 13:24], dmax_mon_2p0K_9[, 13:24],
                               dmax_mon_2p0K_10[, 13:24], dmax_mon_2p0K_11[, 13:24], dmax_mon_2p0K_12[, 13:24],
                               dmax_mon_2p0K_13[, 13:24])
max_dis_mon_2p0K_base <- rbind(dmax_mon_2p0K_1[, 25:36], dmax_mon_2p0K_2[, 25:36], dmax_mon_2p0K_3[, 25:36],
                               dmax_mon_2p0K_4[, 25:36], dmax_mon_2p0K_5[, 25:36], dmax_mon_2p0K_6[, 25:36],
                               dmax_mon_2p0K_7[, 25:36], dmax_mon_2p0K_8[, 25:36], dmax_mon_2p0K_9[, 25:36],
                               dmax_mon_2p0K_10[, 25:36], dmax_mon_2p0K_11[, 25:36], dmax_mon_2p0K_12[, 25:36],
                               dmax_mon_2p0K_13[, 25:36])

max_mel_mon_2p0K_koel <- rbind(fmax_mon_2p0K_1[, 25:36], fmax_mon_2p0K_2[, 25:36], fmax_mon_2p0K_3[, 25:36],
                               fmax_mon_2p0K_4[, 25:36], fmax_mon_2p0K_5[, 25:36], fmax_mon_2p0K_6[, 25:36],
                               fmax_mon_2p0K_7[, 25:36], fmax_mon_2p0K_8[, 25:36], fmax_mon_2p0K_9[, 25:36],
                               fmax_mon_2p0K_10[, 25:36], fmax_mon_2p0K_11[, 25:36], fmax_mon_2p0K_12[, 25:36],
                               fmax_mon_2p0K_13[, 25:36])
max_mel_mon_2p0K_coch <- rbind(fmax_mon_2p0K_1[, 13:24], fmax_mon_2p0K_2[, 13:24], fmax_mon_2p0K_3[, 13:24],
                               fmax_mon_2p0K_4[, 13:24], fmax_mon_2p0K_5[, 13:24], fmax_mon_2p0K_6[, 13:24],
                               fmax_mon_2p0K_7[, 13:24], fmax_mon_2p0K_8[, 13:24], fmax_mon_2p0K_9[, 13:24],
                               fmax_mon_2p0K_10[, 13:24], fmax_mon_2p0K_11[, 13:24], fmax_mon_2p0K_12[, 13:24],
                               fmax_mon_2p0K_13[, 13:24])
max_mel_mon_2p0K_base <- rbind(fmax_mon_2p0K_1[, 1:12], fmax_mon_2p0K_2[, 1:12], fmax_mon_2p0K_3[, 1:12],
                               fmax_mon_2p0K_4[, 1:12], fmax_mon_2p0K_5[, 1:12], fmax_mon_2p0K_6[, 1:12],
                               fmax_mon_2p0K_7[, 1:12], fmax_mon_2p0K_8[, 1:12], fmax_mon_2p0K_9[, 1:12],
                               fmax_mon_2p0K_10[, 1:12], fmax_mon_2p0K_11[, 1:12], fmax_mon_2p0K_12[, 1:12],
                               fmax_mon_2p0K_13[, 1:12])

max_prt_mon_2p0K_koel <- rbind(pmax_mon_2p0K_1[, 25:36], pmax_mon_2p0K_2[, 25:36], pmax_mon_2p0K_3[, 25:36],
                               pmax_mon_2p0K_4[, 25:36], pmax_mon_2p0K_5[, 25:36], pmax_mon_2p0K_6[, 25:36],
                               pmax_mon_2p0K_7[, 25:36], pmax_mon_2p0K_8[, 25:36], pmax_mon_2p0K_9[, 25:36],
                               pmax_mon_2p0K_10[, 25:36], pmax_mon_2p0K_11[, 25:36], pmax_mon_2p0K_12[, 25:36],
                               pmax_mon_2p0K_13[, 25:36])
max_prt_mon_2p0K_coch <- rbind(pmax_mon_2p0K_1[, 13:24], pmax_mon_2p0K_2[, 13:24], pmax_mon_2p0K_3[, 13:24],
                               pmax_mon_2p0K_4[, 13:24], pmax_mon_2p0K_5[, 13:24], pmax_mon_2p0K_6[, 13:24],
                               pmax_mon_2p0K_7[, 13:24], pmax_mon_2p0K_8[, 13:24], pmax_mon_2p0K_9[, 13:24],
                               pmax_mon_2p0K_10[, 13:24], pmax_mon_2p0K_11[, 13:24], pmax_mon_2p0K_12[, 13:24],
                               pmax_mon_2p0K_13[, 13:24])
max_prt_mon_2p0K_base <- rbind(pmax_mon_2p0K_1[, 1:12], pmax_mon_2p0K_2[, 1:12], pmax_mon_2p0K_3[, 1:12],
                               pmax_mon_2p0K_4[, 1:12], pmax_mon_2p0K_5[, 1:12], pmax_mon_2p0K_6[, 1:12],
                               pmax_mon_2p0K_7[, 1:12], pmax_mon_2p0K_8[, 1:12], pmax_mon_2p0K_9[, 1:12],
                               pmax_mon_2p0K_10[, 1:12], pmax_mon_2p0K_11[, 1:12], pmax_mon_2p0K_12[, 1:12],
                               pmax_mon_2p0K_13[, 1:12])

max_prl_mon_2p0K_koel <- rbind(pmax_mon_2p0K_1[, 61:72], pmax_mon_2p0K_2[, 61:72], pmax_mon_2p0K_3[, 61:72],
                               pmax_mon_2p0K_4[, 61:72], pmax_mon_2p0K_5[, 61:72], pmax_mon_2p0K_6[, 61:72],
                               pmax_mon_2p0K_7[, 61:72], pmax_mon_2p0K_8[, 61:72], pmax_mon_2p0K_9[, 61:72],
                               pmax_mon_2p0K_10[, 61:72], pmax_mon_2p0K_11[, 61:72], pmax_mon_2p0K_12[, 61:72],
                               pmax_mon_2p0K_13[, 61:72])
max_prl_mon_2p0K_coch <- rbind(pmax_mon_2p0K_1[, 49:60], pmax_mon_2p0K_2[, 49:60], pmax_mon_2p0K_3[, 49:60],
                               pmax_mon_2p0K_4[, 49:60], pmax_mon_2p0K_5[, 49:60], pmax_mon_2p0K_6[, 49:60],
                               pmax_mon_2p0K_7[, 49:60], pmax_mon_2p0K_8[, 49:60], pmax_mon_2p0K_9[, 49:60],
                               pmax_mon_2p0K_10[, 49:60], pmax_mon_2p0K_11[, 49:60], pmax_mon_2p0K_12[, 49:60],
                               pmax_mon_2p0K_13[, 49:60])
max_prl_mon_2p0K_base <- rbind(pmax_mon_2p0K_1[, 37:48], pmax_mon_2p0K_2[, 37:48], pmax_mon_2p0K_3[, 37:48],
                               pmax_mon_2p0K_4[, 37:48], pmax_mon_2p0K_5[, 37:48], pmax_mon_2p0K_6[, 37:48],
                               pmax_mon_2p0K_7[, 37:48], pmax_mon_2p0K_8[, 37:48], pmax_mon_2p0K_9[, 37:48],
                               pmax_mon_2p0K_10[, 37:48], pmax_mon_2p0K_11[, 37:48], pmax_mon_2p0K_12[, 37:48],
                               pmax_mon_2p0K_13[, 37:48])

#Put together 3.0K warming level
max_dis_mon_3p0K_koel <- rbind(dmax_mon_3p0K_1[, 1:12], dmax_mon_3p0K_2[, 1:12], dmax_mon_3p0K_3[, 1:12],
                               dmax_mon_3p0K_4[, 1:12], dmax_mon_3p0K_5[, 1:12], dmax_mon_3p0K_6[, 1:12],
                               dmax_mon_3p0K_7[, 1:12], dmax_mon_3p0K_8[, 1:12])
max_dis_mon_3p0K_coch <- rbind(dmax_mon_3p0K_1[, 13:24], dmax_mon_3p0K_2[, 13:24], dmax_mon_3p0K_3[, 13:24],
                               dmax_mon_3p0K_4[, 13:24], dmax_mon_3p0K_5[, 13:24], dmax_mon_3p0K_6[, 13:24],
                               dmax_mon_3p0K_7[, 13:24], dmax_mon_3p0K_8[, 13:24])
max_dis_mon_3p0K_base <- rbind(dmax_mon_3p0K_1[, 25:36], dmax_mon_3p0K_2[, 25:36], dmax_mon_3p0K_3[, 25:36],
                               dmax_mon_3p0K_4[, 25:36], dmax_mon_3p0K_5[, 25:36], dmax_mon_3p0K_6[, 25:36],
                               dmax_mon_3p0K_7[, 25:36], dmax_mon_3p0K_8[, 25:36])

max_mel_mon_3p0K_koel <- rbind(fmax_mon_3p0K_1[, 25:36], fmax_mon_3p0K_2[, 25:36], fmax_mon_3p0K_3[, 25:36],
                               fmax_mon_3p0K_4[, 25:36], fmax_mon_3p0K_5[, 25:36], fmax_mon_3p0K_6[, 25:36],
                               fmax_mon_3p0K_7[, 25:36], fmax_mon_3p0K_8[, 25:36])
max_mel_mon_3p0K_coch <- rbind(fmax_mon_3p0K_1[, 13:24], fmax_mon_3p0K_2[, 13:24], fmax_mon_3p0K_3[, 13:24],
                               fmax_mon_3p0K_4[, 13:24], fmax_mon_3p0K_5[, 13:24], fmax_mon_3p0K_6[, 13:24],
                               fmax_mon_3p0K_7[, 13:24], fmax_mon_3p0K_8[, 13:24])
max_mel_mon_3p0K_base <- rbind(fmax_mon_3p0K_1[, 1:12], fmax_mon_3p0K_2[, 1:12], fmax_mon_3p0K_3[, 1:12],
                               fmax_mon_3p0K_4[, 1:12], fmax_mon_3p0K_5[, 1:12], fmax_mon_3p0K_6[, 1:12],
                               fmax_mon_3p0K_7[, 1:12], fmax_mon_3p0K_8[, 1:12])

max_prt_mon_3p0K_koel <- rbind(pmax_mon_3p0K_1[, 25:36], pmax_mon_3p0K_2[, 25:36], pmax_mon_3p0K_3[, 25:36],
                               pmax_mon_3p0K_4[, 25:36], pmax_mon_3p0K_5[, 25:36], pmax_mon_3p0K_6[, 25:36],
                               pmax_mon_3p0K_7[, 25:36], pmax_mon_3p0K_8[, 25:36])
max_prt_mon_3p0K_coch <- rbind(pmax_mon_3p0K_1[, 13:24], pmax_mon_3p0K_2[, 13:24], pmax_mon_3p0K_3[, 13:24],
                               pmax_mon_3p0K_4[, 13:24], pmax_mon_3p0K_5[, 13:24], pmax_mon_3p0K_6[, 13:24],
                               pmax_mon_3p0K_7[, 13:24], pmax_mon_3p0K_8[, 13:24])
max_prt_mon_3p0K_base <- rbind(pmax_mon_3p0K_1[, 1:12], pmax_mon_3p0K_2[, 1:12], pmax_mon_3p0K_3[, 1:12],
                               pmax_mon_3p0K_4[, 1:12], pmax_mon_3p0K_5[, 1:12], pmax_mon_3p0K_6[, 1:12],
                               pmax_mon_3p0K_7[, 1:12], pmax_mon_3p0K_8[, 1:12])

max_prl_mon_3p0K_koel <- rbind(pmax_mon_3p0K_1[, 61:72], pmax_mon_3p0K_2[, 61:72], pmax_mon_3p0K_3[, 61:72],
                               pmax_mon_3p0K_4[, 61:72], pmax_mon_3p0K_5[, 61:72], pmax_mon_3p0K_6[, 61:72],
                               pmax_mon_3p0K_7[, 61:72], pmax_mon_3p0K_8[, 61:72])
max_prl_mon_3p0K_coch <- rbind(pmax_mon_3p0K_1[, 49:60], pmax_mon_3p0K_2[, 49:60], pmax_mon_3p0K_3[, 49:60],
                               pmax_mon_3p0K_4[, 49:60], pmax_mon_3p0K_5[, 49:60], pmax_mon_3p0K_6[, 49:60],
                               pmax_mon_3p0K_7[, 49:60], pmax_mon_3p0K_8[, 49:60])
max_prl_mon_3p0K_base <- rbind(pmax_mon_3p0K_1[, 37:48], pmax_mon_3p0K_2[, 37:48], pmax_mon_3p0K_3[, 37:48],
                               pmax_mon_3p0K_4[, 37:48], pmax_mon_3p0K_5[, 37:48], pmax_mon_3p0K_6[, 37:48],
                               pmax_mon_3p0K_7[, 37:48], pmax_mon_3p0K_8[, 37:48])

plot_month_max <- function(max_hist_mon, max_1p5K_mon, max_2p0K_mon, max_3p0K_mon,
                           ylims = c(0, 100), ylab = "", cex_axs = 1.6, cex_lab = 1.5,
                           main = "", cex_main = 1.8, do_legend = F, legend_pos = "topright"){
  
  col_hist <- "steelblue4"
  col_1p5K <- "grey25"
  col_2p0K <- "orange3"
  col_3p0K <- "darkred"
  
  x_labs <- c("O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A", "S")
  x_ats <- seq(2.5, 60, 5)
  
  plot(1:10, 1:10, type = "n", ylim = ylims, xlim = c(1, 60), axes = F,
       xlab = "", ylab = "")
  grid(nx = 0, ny = 6, col = "grey55", lwd = 0.5)
  boxplot(max_hist_mon, at = seq(1, 60, 5), add = T, col = col_hist,
          outline = F, notch = T, axes = F, whisklty = 0, staplelty = 0)
  boxplot(max_1p5K_mon, at = seq(2, 60, 5), add = T, col = col_1p5K,
          outline = F, notch = T, axes = F, whisklty = 0, staplelty = 0)
  boxplot(max_2p0K_mon, at = seq(3, 60, 5), add = T, col = col_2p0K,
          outline = F, notch = T, axes = F, whisklty = 0, staplelty = 0)
  boxplot(max_3p0K_mon, at = seq(4, 60, 5), add = T, col = col_3p0K,
          outline = F, notch = T, axes = F, whisklty = 0, staplelty = 0)
  axis(2, mgp=c(3, 0.35, 0), tck = -0.015, cex.axis = cex_axs)
  axis(1, labels = x_labs, at = x_ats, mgp=c(3, 0.42, 0), tck = -0.0, cex.axis = cex_axs)
  axis(1, labels = F, at = seq(5, 55, 5), tck = -0.045)
  mtext(main, side = 3, line = 0.25, adj = 0, cex = cex_main)
  mtext(ylab, side = 2, line = 2.2, cex = cex_lab)
  if(do_legend){
    legend(legend_pos, c("Hist.", "1.5K", "2.0K", "3.0K"), pch = 19, 
           col = c(col_hist, col_1p5K, col_2p0K, col_3p0K), cex = 1.5,
           box.lwd = 0.0, box.col = "black", bg = "white", ncol = 2)
  }
  box()
  
}

#Plot: Monthly maxima discharge
pdf(paste0(bas_dir,"res_figs/mon_dis_fut.pdf"), width = 8, height = 8)

par(family = "serif")
par(mar = c(1.5, 4.8, 3.0, 1.0))

layout(matrix(c(1, 2, 3),
              3, 1, byrow = T), widths=c(), heights=c())
# layout.show(n = 8)

plot_month_max(max_dis_mon_hist_base, max_dis_mon_1p5K_base, 
               max_dis_mon_2p0K_base, max_dis_mon_3p0K_base, ylims = c(800, 2800),
               main = "a) Basel", ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")),
               do_legend = T)

plot_month_max(max_dis_mon_hist_coch, max_dis_mon_1p5K_coch, 
               max_dis_mon_2p0K_coch, max_dis_mon_3p0K_coch, ylims = c(80, 1600),
               main = "b) Cochem", ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))

plot_month_max(max_dis_mon_hist_koel, max_dis_mon_1p5K_koel, 
               max_dis_mon_2p0K_koel, max_dis_mon_3p0K_koel, ylims = c(1300, 6600),
               main = "c) Cologne", ylab = expression(paste("Discharge [m"^"3", "s"^"-1","]")))

dev.off()




#Plot: Monthly maxima discharge
pdf(paste0(bas_dir,"res_figs/mon_add_fut.pdf"), width = 16, height = 8)

par(family = "serif")
par(mar = c(1.5, 3.0, 1.5, 1.0))

layout(matrix(c(7, 7, 7,
                8, 1, 2,
                8, 3, 4,
                8, 5, 6),
              4, 3, byrow = T), widths=c(0.05, 1, 1), heights=c(0.1, 1, 1, 1))
# layout.show(n = 8)
my_cex_axs = 1.8

plot_month_max(max_mel_mon_hist_base, max_mel_mon_1p5K_base, 
               max_mel_mon_2p0K_base, max_mel_mon_3p0K_base, ylims = c(0, 45),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5, do_legend = T)

plot_month_max(max_mel_mon_hist_coch, max_mel_mon_1p5K_coch, 
               max_mel_mon_2p0K_coch, max_mel_mon_3p0K_coch, ylims = c(0, 45),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5)

plot_month_max(max_prl_mon_hist_base, max_prl_mon_1p5K_base, 
               max_prl_mon_2p0K_base, max_prl_mon_3p0K_base, ylims = c(5, 30),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5)

plot_month_max(max_prl_mon_hist_coch, max_prl_mon_1p5K_coch, 
               max_prl_mon_2p0K_coch, max_prl_mon_3p0K_coch, ylims = c(5, 18),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5)

plot_month_max(max_prt_mon_hist_base, max_prt_mon_1p5K_base, 
               max_prt_mon_2p0K_base, max_prt_mon_3p0K_base, ylims = c(5, 30),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5)

plot_month_max(max_prt_mon_hist_coch, max_prt_mon_1p5K_coch, 
               max_prt_mon_2p0K_coch, max_prt_mon_3p0K_coch, ylims = c(5, 18),
               main = "", ylab = "", cex_axs = my_cex_axs, cex_lab = 1.5)

#Gauging station
cex_header <- 1.7
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("a) Basel",
      side = 3, line = -2.55, cex = cex_header+0.2, adj = 0.25)
mtext("b) Cochem",
      side = 3, line = -2.55, cex = cex_header+0.2, adj = 0.80)

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Snowmelt [mm]",  side = 2, line = -2.5, cex = cex_header, adj = 0.90, outer = T)
mtext("2. Prec. liq. [mm]",  side = 2, line = -2.5, cex = cex_header, adj = 0.45, outer = T)
mtext("3. Prec. tot. [mm]",  side = 2, line = -2.5, cex = cex_header, adj = 0.05, outer = T)

dev.off()





#event_analy----

max_na(dmax_1p5K_11[, 1])

which(dmax_1p5K_11[, 1] == max_na(dmax_1p5K_11[, 1]))


#Function to extract time series from nc-files
snow_from_nc <- function(gcm_model, delta_t, rcp, ncores = 5){
  
  #select nc_file
  nc_path_sel <- nc_file_paths[which(grepl(gcm_model, nc_file_paths) &
                                       grepl(rcp, nc_file_paths))]
  nc_file_sel <- nc_open(paste0(nc_path_sel, "output/mHM_Fluxes_States.nc"))
  
  #path to precipitation input
  # nc_path_sel_prec <- gsub("output", "input/meteo", nc_path_sel)
  # pr_file_sel <- nc_open(paste0(nc_path_sel, "input/pre.nc"))
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  
  #if simulation time frame does not entirely cover warming period
  if(date_sel[1] > date[1]){
    sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
    count_date <- length(date_sel)
  }else{
    sta_date_ind <- 1
    count_date <- length(date_sel) - which(format(date_sel) == date[1]) + 1
  }
  
  #snowpack
  snow_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
  
  n_cores <- ncores #number of cores used for parallel computing
  
  #Make cluster for parallel computing
  my_clust <- makeCluster(n_cores)
  clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
  clusterExport(my_clust, c("index_col_base", "index_row_base", "index_col_coch", "index_row_coch",
                            "index_col_koel", "index_row_koel", "elevs_base", "elevs_coch", "elevs_koel"))
  registerDoParallel(my_clust)
  
  sno_base <- foreach(i = 1:length(index_col_base), .combine = 'cbind') %dopar% {
    snow_cube[index_col_base[i], index_row_base[i], ]
  }
  sno_coch <- foreach(i = 1:length(index_col_coch), .combine = 'cbind') %dopar% {
    snow_cube[index_col_coch[i], index_row_coch[i], ]
  }
  
  stopCluster(my_clust)
  
  #Values on basin scale
  base_sd_mea <- apply(sno_base, 1, mea_na)
  coch_sd_mea <- apply(sno_coch, 1, mea_na)
  
  date_out <- date[which(date %in% date_sel)]
  
  snow_out <- data.frame(date = date_out,
                         base_sd_mea = base_sd_mea,
                         coch_sd_mea = coch_sd_mea)  
  
  return(snow_out)
  
}

snow_1p5K_11 <- snow_from_nc("HadGEM2-ES", "1p5", "8p5")
snow_3p0K_1  <- snow_from_nc("HadGEM2-ES", "3p0", "6p0")

disc_sel <- disc_1p5K_11
prec_sel <- prec_1p5K_11
snow_sel <- snow_1p5K_11

ind_sel <- 1370:1550

# disc_sel$date[ind_sel]

plot(disc_sel$koel[ind_sel], type = "l", ylim = c(0, 19000))
lines(disc_sel$base[ind_sel], type = "l", col = "red")
lines(disc_sel$coch[ind_sel], type = "l", col = "blue")

plot(disc_sel$base[ind_sel], type = "l", col = "black")
par(new = T)
plot(prec_sel$base_pre_liq[ind_sel], type = "h", col = "blue2")
par(new = T)
plot(snow_sel$base_sd_mea[ind_sel], type = "l", col = "orange2")

plot(disc_sel$coch[ind_sel], type = "l", col = "black")
par(new = T)
plot(prec_sel$coch_pre_liq[ind_sel], type = "h", col = "blue2")
par(new = T)
plot(snow_sel$coch_sd_mea[ind_sel], type = "l", col = "orange2")



#parde----

#Parde-Coefficients
grdc_data_base <- read_grdc(paste0(grdc_dir, "6935051_Q_Day.Cmd.txt"))
grdc_data_koel <- read_grdc(paste0(grdc_dir, "6335060_Q_Day.Cmd.txt"))
grdc_data_coch <- read_grdc(paste0(grdc_dir, "6336050_Q_Day.Cmd.txt"))

#Order data by day 
data_day_coch <- ord_day(data_in = grdc_data_coch$value,
                         date = grdc_data_coch$date,
                         start_y = 1917,
                         end_y = 2016,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)

data_day_base <- ord_day(data_in = grdc_data_base$value,
                         date = grdc_data_base$date,
                         start_y = 1917,
                         end_y = 2016,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)

data_day_koel <- ord_day(data_in = grdc_data_koel$value,
                         date = grdc_data_koel$date,
                         start_y = 1917,
                         end_y = 2016,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)

yea_mea_base <- apply(data_day_base, 2, mea_na)
yea_mea_koel <- apply(data_day_koel, 2, mea_na)
yea_mea_coch <- apply(data_day_coch, 2, mea_na)

month_ind <- c(rep(1, 31), rep(2, 30), rep(3, 31), rep(4, 31), rep(5, 28), 
               rep(6, 31), rep(7, 30), rep(8, 31), rep(9, 30), rep(10, 31), 
               rep(11, 31), rep(12, 30))#hydrological year: 1 = October

mon_mea_base <- aggregate(yea_mea_base, by = list(months = month_ind), FUN = mea_na)
mon_mea_koel <- aggregate(yea_mea_koel, by = list(months = month_ind), FUN = mea_na)
mon_mea_coch <- aggregate(yea_mea_coch, by = list(months = month_ind), FUN = mea_na)

parde_ind_base <- mon_mea_base$x / mea_na(yea_mea_base)
parde_ind_koel <- mon_mea_koel$x / mea_na(yea_mea_koel)
parde_ind_coch <- mon_mea_coch$x / mea_na(yea_mea_coch)

pdf(paste0(bas_dir,"res_figs/parde_dis.pdf"), width = 8, height = 3)

ylims <- range(c(parde_ind_coch, parde_ind_base, parde_ind_koel))

x_axis_lab <- 1:12
x_axis_tic <- (1:13)-0.5
col_coch <- "steelblue4" 
col_base <- "darkred"
col_koel <- "black"

par(family = "serif")
par(mar = c(1.8, 2.6, 1.5, 0.2))

plot(parde_ind_coch, type = "n", axes = F, ylim = ylims, ylab = "", xlab = "")
abline(v = c(x_axis_tic), lty = "dashed", col = "grey55", lwd = 0.8)
lines(parde_ind_base, col = col_base, lwd = 2.2)
lines(parde_ind_koel, col = col_koel, lwd = 2.2)
lines(parde_ind_coch, col = col_coch, lwd = 2.2)
points(parde_ind_coch, col = col_coch, pch = 19, cex = 0.8)
points(parde_ind_base, col = col_base, pch = 19, cex = 0.8)
points(parde_ind_koel, col = col_koel, pch = 19, cex = 0.8)
axis(2, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.35, 0), cex.axis = 1.2)#plot labels
mtext("Pardé-coefficient", side = 2, line = 1.3, cex = 1.4, adj = 0.5)
mtext("Runoff seasonality", side = 3, line = 0.3, cex = 1.6, adj = 0.0)
legend("topright", c("Cochem", "Cologne", "Basel"), pch = 19, ncol = 3, bg = "white",
       col = c(col_coch, col_koel, col_base))
box()

dev.off()


#regi_illus----

pdf(paste0(bas_dir,"res_figs/regi_illus.pdf"), width = 8, height = 3.5)

#Positions ticks and labels for x-axis
x_axis_lab <- c(15,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(15,46,74,105,135,166,196,227,258,288,319,349,379)-14

par(oma=c(0,0,0,0))
par(mar = c(1.5, 0.5, 2.0, 0.5))
par(family = "serif")

col_niv <- "darkred" 
col_plu <- "steelblue4"

plot(1:200, 1:200, type="n", ylim = c(0, 0.9), xlim = c(1, 365), axes=F, xaxs="i",yaxs="i")

axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.5, 0), cex.axis = 1.5)
#pluvial
x   <- seq(1, 210, length=1000)
y   <- dnorm(x, mean = 109, sd = 34)
lines(x, y*69, type="l", lty="dashed", col = col_plu, lwd=2.3)
lines(x, y*55, type="l", lty="solid", col = col_plu, lwd=2.3)
#nival
x   <- seq(92+75,230+75,length=1000)
y   <- dnorm(x, mean=155+75, sd=22)
lines(x, y*35, type="l", lty="dashed", col = col_niv, lwd=2.3)
x   <- seq(88+118,225+118,length=1000)
y   <- dnorm(x, mean=150+118, sd=22)
lines(x, y*35, type="l", lty="solid", col = col_niv, lwd=2.3)
#overlap
x   <- seq(85+67,155+71,length=1000)
y   <- dnorm(x, mean=122+63, sd=13)
polygon(x, y*1.9, col="black")  
#arrows
shape::Arrows(x0=210+58, x1=168+62, y0=0.68, y1=0.68, lwd=2, col=col_niv)#nival
shape::Arrows(x0=109, x1=109, y0=0.67, y1=0.76, lwd=2, col=col_plu)#pluvial
#Question marks
text(140, 0.76, "?", cex=2, col=col_plu)
text(188+61, 0.76, "?", cex=2, col=col_niv)
#Text
text(42, 0.62, "pluvial", cex=2, col=col_plu)
text(305, 0.62, "nival", cex=2, col=col_niv)
text(140, 0.08, "overlap?", cex=2, col="black")
mtext("Idealized seasonal distribution of extreme discharge", side = 3, line = 0.3, cex = 1.6, adj = 0.0)
#Legend
legend("topright", legend = c("present", "future"), lty = c("solid","dashed"),
       bty = "n", cex=1.4)
# mtext("Idealized seasonal distribution extreme discharges", side = 3, cex=1.5, line = 0.2)

box()

dev.off()
