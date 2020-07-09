###

#mHM simulations Part II: mHM simulations future scenarios RCP2.6, RCP6.0 and RCP8.5
#Erwin Rottler, UniversityI of Potsdam, Summer 2020

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales, ncdf4,
               ncdf4.helpers, sp, raster, viridis, meltimr, POT, readr, hydroGOF,
               CoinCalc, seas, tictoc)

#set directories
bas_dir <- "U:/rhine_fut/R/"
run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

#load functions
source(paste0(bas_dir, "mhm_rhine/functs.R"))

stopCluster(my_clust)

n_cores <- 15 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, raster))
registerDoParallel(my_clust)

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#read table warming levels
warming_periods <- read.table("D:/nrc_user/rottler/mhm_run/6435060/input/meteo/periods_warming.txt",
                              header = T, stringsAsFactors = F)


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



#nc_files----


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
  
  return(disc_sim)
  
}

#historical
disc_hist_1 <- disc_from_nc("GFDL-ESM2M",     "historical", "historical")
disc_hist_2 <- disc_from_nc("HadGEM2-ES",     "historical", "historical")
disc_hist_3 <- disc_from_nc("IPSL-CM5A-LR",   "historical", "historical")
disc_hist_4 <- disc_from_nc("MIROC-ESM-CHEM", "historical", "historical")
disc_hist_5 <- disc_from_nc("NorESM1-M",      "historical", "historical")

#1.5K warming level
disc_1p5K_1  <- disc_from_nc("HadGEM2-ES",     "1p5", "2p6")
disc_1p5K_2  <- disc_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
disc_1p5K_3  <- disc_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
disc_1p5K_3  <- disc_from_nc("NorESM1-M",      "1p5", "2p6")
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

#2K warming level
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

#3K warming level
disc_3p0K_1 <- disc_from_nc("HadGEM2-ES",     "3p0", "6p0")
disc_3p0K_2 <- disc_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
disc_3p0K_3 <- disc_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
disc_3p0K_4 <- disc_from_nc("GFDL-ESM2M",     "3p0", "8p5")
disc_3p0K_5 <- disc_from_nc("HadGEM2-ES",     "3p0", "8p5")
disc_3p0K_6 <- disc_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
disc_3p0K_7 <- disc_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
disc_3p0K_8 <- disc_from_nc("NorESM1-M",      "3p0", "8p5")


#get_fluxes----

#Get basin .shp
basin_coch_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/cochem_catch.shp")
basin_base_raw <- rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/eu_dem/processed/basins/basel_catch.shp")
basin_coch <- spTransform(basin_coch_raw, CRS = crswgs84)
basin_base <- spTransform(basin_base_raw, CRS = crswgs84)

#get cells in basin
lon <- ncdf4::ncvar_get(nc_cm1_disc_hist, varid = "lon")
lat <- ncdf4::ncvar_get(nc_cm1_disc_hist, varid = "lat")

#spatial grid points from lat/lon
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)

#grid points inside watersheds
inside_base <- !is.na(sp::over(grid_points_cube_84, as(basin_base, "SpatialPolygons")))
inside_coch <- !is.na(sp::over(grid_points_cube_84, as(basin_coch, "SpatialPolygons")))
grid_points_base <- grid_points_cube_84[which(inside_base == T)]
grid_points_coch <- grid_points_cube_84[which(inside_coch == T)]

#Select cells for Basel/Cochem watershed
lat_in_base <- grid_points_base@coords[, 2]
lat_in_coch <- grid_points_coch@coords[, 2]

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

#paths to mHM nc-files with fluxes and states
nc_file_flux_paths <- 
  c(paste0(run_dir, "output/GFDL-ESM2M/historical/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/historical/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/historical/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/historical/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/historical/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp2p6/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp2p6/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp2p6/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp2p6/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp2p6/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp6p0/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp6p0/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp6p0/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp6p0/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp6p0/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/GFDL-ESM2M/rcp8p5/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/HadGEM2-ES/rcp8p5/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/IPSL-CM5A-LR/rcp8p5/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/MIROC-ESM-CHEM/rcp8p5/output/mHM_Fluxes_States.nc"),
    paste0(run_dir, "output/NorESM1-M/rcp8p5/output/mHM_Fluxes_States.nc"))



flux_from_nc <- function(gcm_model, delta_t, rcp){
  
  #select nc_file
  nc_path_sel <- nc_file_flux_paths[which(grepl(gcm_model, nc_file_flux_paths) &
                                            grepl(rcp, nc_file_flux_paths))]
  nc_file_sel <- nc_open(nc_path_sel)
  
  #get warming period
  wp_years <- get_warming_period(gcm_model, delta_t, rcp)
  
  date_sel <- seq(as.Date(paste0(wp_years[1], "-01-01"), format = "%Y-%m-%d"), 
                  as.Date(paste0(wp_years[2], "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  #date from nc-file
  date <- as.Date(as.character(nc.get.time.series(nc_file_sel, time.dim.name = "time")))
  sta_date_ind <- which(format(date) == paste0(wp_years[1], "-01-01"))
  count_date <- length(date_sel)
  
  #snowpack
  snow_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                         count = c(nrow(lon), ncol(lon), count_date), varid = "snowpack")
  #effective precipitation
  pef_cube <- ncvar_get(nc_file_sel, start = c(1, 1, sta_date_ind), 
                        count = c(nrow(lon), ncol(lon), count_date), varid = "preEffect")
  
  #get simulation results for Basel watershed
  for (i in 1:length(index_col_base)) {
    
    sno_sing <- snow_cube[cube_index_col_base[i], cube_index_row_base[i], ]
    epn_sing <- pef_cube [cube_index_col_base[i], cube_index_row_base[i], ]
    
    if(i == 1){
      
      sno_base <- sno_sing
      epn_base <- epn_sing
      
    }else{
      
      sno_base <- cbind(sno_base, sno_sing)
      epn_base <- cbind(epn_base, epn_sing)
      
    }
    
  }
  
  #get simulation results for Cochem watershed
  for (i in 1:length(index_col_coch)) {
    
    sno_sing <- snow_cube[cube_index_col_coch[i], cube_index_row_coch[i], ]
    epn_sing <- pef_cube [cube_index_col_coch[i], cube_index_row_coch[i], ]
    
    if(i == 1){
      
      sno_coch <- sno_sing
      epn_coch <- epn_sing
      
    }else{
      
      sno_coch <- cbind(sno_coch, sno_sing)
      epn_coch <- cbind(epn_coch, epn_sing)
      
    }
    
  }
  
  #Values on basin scale
  base_sd_mea <- c(NA, apply(sno_base, 1, mea_na))
  base_sd_mea_dif <- c(NA, diff(base_sd_mea))
  base_sd_mea_dif[which(base_sd_mea_dif > 0)] <- NA
  base_me_mea <- base_sd_mea_dif * -1 #melt positive values
  base_ep_mea <- c(NA, apply(epn_base, 1, mea_na))
  
  coch_sd_mea <- c(NA, apply(sno_coch, 1, mea_na))
  coch_sd_mea_dif <- c(NA, diff(coch_sd_mea))
  coch_sd_mea_dif[which(coch_sd_mea_dif > 0)] <- NA
  coch_me_mea <- coch_sd_mea_dif * -1 #melt positive values
  coch_ep_mea <- c(NA, apply(epn_coch, 1, mea_na)) #fluxes/states only start from 02.01.1954
  
  base_melt_ma <- rollapply(data = base_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
  coch_melt_ma <- rollapply(data = coch_me_mea, width = 14,
                            FUN = sum_na, align = "right", fill = NA)
  
  #Moving window liquid precipitation
  base_lpre <- base_ep_mea - base_me_mea
  coch_lpre <- coch_ep_mea - coch_me_mea
  base_lpre_ma <- rollapply(data = base_lpre, width = 5,
                            FUN = sum_na, align = "right", fill = NA)
  coch_lpre_ma <- rollapply(data = coch_lpre, width = 5,
                            FUN = sum_na, align = "right", fill = NA)
  
  #fraction snowmelt contribution
  base_frac_me <- base_melt_ma / (base_melt_ma + base_lpre_ma)
  coch_frac_me <- coch_melt_ma / (coch_melt_ma + coch_lpre_ma)
  
  flux_out <- data.frame(base_melt_ma = base_melt_ma,
                         coch_melt_ma = coch_melt_ma,
                         base_lpre_ma = base_lpre_ma,
                         coch_lpre_ma = coch_lpre_ma,
                         base_frac_me = base_frac_me,
                         coch_frac_me = coch_frac_me)  
  
  return(flux_out)
  
}

#historical
tic()
flux_hist_1 <- flux_from_nc("GFDL-ESM2M",     "historical", "historical")
flux_hist_2 <- flux_from_nc("HadGEM2-ES",     "historical", "historical")
flux_hist_3 <- flux_from_nc("IPSL-CM5A-LR",   "historical", "historical")
flux_hist_4 <- flux_from_nc("MIROC-ESM-CHEM", "historical", "historical")
flux_hist_5 <- flux_from_nc("NorESM1-M",      "historical", "historical")
toc()

#1.5K warming level
flux_1p5K_1  <- flux_from_nc("HadGEM2-ES",     "1p5", "2p6")
flux_1p5K_2  <- flux_from_nc("IPSL-CM5A-LR",   "1p5", "2p6")
flux_1p5K_3  <- flux_from_nc("MIROC-ESM-CHEM", "1p5", "2p6")
flux_1p5K_3  <- flux_from_nc("NorESM1-M",      "1p5", "2p6")
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

#2K warming level
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

#3K warming level
flux_3p0K_1 <- flux_from_nc("HadGEM2-ES",     "3p0", "6p0")
flux_3p0K_2 <- flux_from_nc("IPSL-CM5A-LR",   "3p0", "6p0")
flux_3p0K_3 <- flux_from_nc("MIROC-ESM-CHEM", "3p0", "6p0")
flux_3p0K_4 <- flux_from_nc("GFDL-ESM2M",     "3p0", "8p5")
flux_3p0K_5 <- flux_from_nc("HadGEM2-ES",     "3p0", "8p5")
flux_3p0K_6 <- flux_from_nc("IPSL-CM5A-LR",   "3p0", "8p5")
flux_3p0K_7 <- flux_from_nc("MIROC-ESM-CHEM", "3p0", "8p5")
flux_3p0K_8 <- flux_from_nc("NorESM1-M",      "3p0", "8p5")

#ann_max----

#get annual maxima discharge characteristics
f_max_mag_doy <- function(data_in, start_y, end_y, break_day = 274){
  
  date <- seq(as.Date(paste0(start_y, "-01-01"), format = "%Y-%m-%d"), 
              as.Date(paste0(end_y, "-12-31"), format = "%Y-%m-%d"), by = "day")
  
  max_out <- NULL
  for(i in 1:ncol(data_in)){
    
    data_day <- meltimr::ord_day(data_in = data_in[, i],
                                 date = date,
                                 start_y = start_y,
                                 end_y = end_y,
                                 break_day = break_day)
    
    data_max <- apply(data_day, 1, max_na)
    
    max_tim <- function(data_in){
      
      doy_max <- which(data_in == max_na(data_in))
      
      return(doy_max)
      
    }
    
    doy_max <- apply(data_day, 1, max_tim)
    
    max_out <- cbind(max_out, data_max, doy_max)
    
  }
  
  return(max_out)
  
}

#Annual maxima discharge: Historical
dmax_hist_1 <- f_max_mag_doy(disc_hist_1, 1971, 2000)
dmax_hist_2 <- f_max_mag_doy(disc_hist_2, 1971, 2000)
dmax_hist_3 <- f_max_mag_doy(disc_hist_3, 1971, 2000)
dmax_hist_4 <- f_max_mag_doy(disc_hist_4, 1971, 2000)
dmax_hist_5 <- f_max_mag_doy(disc_hist_5, 1971, 2000)
fmax_hist_1 <- f_max_mag_doy(flux_hist_1, 1971, 2000)
fmax_hist_2 <- f_max_mag_doy(flux_hist_2, 1971, 2000)
fmax_hist_3 <- f_max_mag_doy(flux_hist_3, 1971, 2000)
fmax_hist_4 <- f_max_mag_doy(flux_hist_4, 1971, 2000)
fmax_hist_5 <- f_max_mag_doy(flux_hist_5, 1971, 2000)

#Annual maxima discharge: 3K warming
dmax_3K_1 <- f_max_mag_doy(disc_3K_1, 2067, 2096)
dmax_3K_2 <- f_max_mag_doy(disc_3K_2, 2035, 2064)
dmax_3K_3 <- f_max_mag_doy(disc_3K_3, 2038, 2067)
dmax_3K_4 <- f_max_mag_doy(disc_3K_4, 2037, 2066)
dmax_3K_5 <- f_max_mag_doy(disc_3K_5, 2057, 2086)
dmax_3K_6 <- f_max_mag_doy(disc_3K_6, 2056, 2085)
dmax_3K_7 <- f_max_mag_doy(disc_3K_7, 2066, 2095)
dmax_3K_8 <- f_max_mag_doy(disc_3K_8, 2055, 2084)
fmax_3K_1 <- f_max_mag_doy(flux_3K_1, 2067, 2096)
fmax_3K_2 <- f_max_mag_doy(flux_3K_2, 2035, 2064)
fmax_3K_3 <- f_max_mag_doy(flux_3K_3, 2038, 2067)
fmax_3K_4 <- f_max_mag_doy(flux_3K_4, 2037, 2066)
fmax_3K_5 <- f_max_mag_doy(flux_3K_5, 2057, 2086)
fmax_3K_6 <- f_max_mag_doy(flux_3K_6, 2056, 2085)
fmax_3K_7 <- f_max_mag_doy(flux_3K_7, 2066, 2095)
fmax_3K_8 <- f_max_mag_doy(flux_3K_8, 2055, 2084)

#Annual maxima discharge: 2K warming
dmax_2K_1 <- f_max_mag_doy(disc_2K_1, 2038, 2067)
dmax_2K_2 <- f_max_mag_doy(disc_2K_2, 2016, 2045)
dmax_2K_3 <- f_max_mag_doy(disc_2K_3, 2018, 2047)
dmax_2K_4 <- f_max_mag_doy(disc_2K_4, 2017, 2046)
dmax_2K_5 <- f_max_mag_doy(disc_2K_5, 2031, 2060)
dmax_2K_6 <- f_max_mag_doy(disc_2K_6, 2060, 2089)
dmax_2K_7 <- f_max_mag_doy(disc_2K_7, 2026, 2055)
dmax_2K_8 <- f_max_mag_doy(disc_2K_8, 2028, 2057)
dmax_2K_9 <- f_max_mag_doy(disc_2K_9, 2028, 2057)
dmax_2K_10 <- f_max_mag_doy(disc_2K_10, 2054, 2083)
dmax_2K_11 <- f_max_mag_doy(disc_2K_11, 2029, 2058)
dmax_2K_12 <- f_max_mag_doy(disc_2K_12, 2060, 2089)
dmax_2K_13 <- f_max_mag_doy(disc_2K_13, 2023, 2052)

#Annual maxima discharge: 1.5K warming
dmax_1p5K_1 <- f_max_mag(disc_1p5K_1, 2021, 2050)
dmax_1p5K_2 <- f_max_mag(disc_1p5K_2, 2005, 2033) #simulation only from 2005 on
dmax_1p5K_3 <- f_max_mag(disc_1p5K_3, 2006, 2035)
dmax_1p5K_4 <- f_max_mag(disc_1p5K_4, 2006, 2035)
dmax_1p5K_5 <- f_max_mag(disc_1p5K_5, 2016, 2045)
dmax_1p5K_6 <- f_max_mag(disc_1p5K_6, 2040, 2069)
dmax_1p5K_7 <- f_max_mag(disc_1p5K_7, 2011, 2040)
dmax_1p5K_8 <- f_max_mag(disc_1p5K_8, 2009, 2038)
dmax_1p5K_9 <- f_max_mag(disc_1p5K_9, 2012, 2041)
dmax_1p5K_10 <- f_max_mag(disc_1p5K_10, 2031, 2060)
dmax_1p5K_11 <- f_max_mag(disc_1p5K_11, 2007, 2036)
dmax_1p5K_12 <- f_max_mag(disc_1p5K_12, 2008, 2037)
dmax_1p5K_13 <- f_max_mag(disc_1p5K_13, 2006, 2035)
dmax_1p5K_14 <- f_max_mag(disc_1p5K_14, 2047, 2076)



ord_day(flux_hist_1)

max_fra_hist_1



col_sel <- 2
color_3K <- "red3"
color_2K <- "orange2"
color_1p5K <- "grey25"
color_hist <- "blue3"


plot(max_doy_3K_1[, col_sel], type = "n", max_mag_3K_1[, col_sel], xlim = c(1, 365), ylim = c(0, 5000))

points(max_doy_3K_1[, col_sel], max_mag_3K_1[, col_sel], col = color_3K)
points(max_doy_3K_2[, col_sel], max_mag_3K_2[, col_sel], col = color_3K)
points(max_doy_3K_3[, col_sel], max_mag_3K_3[, col_sel], col = color_3K)
points(max_doy_3K_4[, col_sel], max_mag_3K_4[, col_sel], col = color_3K)
points(max_doy_3K_5[, col_sel], max_mag_3K_5[, col_sel], col = color_3K)
points(max_doy_3K_6[, col_sel], max_mag_3K_6[, col_sel], col = color_3K)
points(max_doy_3K_7[, col_sel], max_mag_3K_7[, col_sel], col = color_3K)
points(max_doy_3K_8[, col_sel], max_mag_3K_8[, col_sel], col = color_3K)

points(max_doy_2K_1[, col_sel], max_mag_2K_1[, col_sel], col = color_2K)
points(max_doy_2K_2[, col_sel], max_mag_2K_2[, col_sel], col = color_2K)
points(max_doy_2K_3[, col_sel], max_mag_2K_3[, col_sel], col = color_2K)
points(max_doy_2K_4[, col_sel], max_mag_2K_4[, col_sel], col = color_2K)
points(max_doy_2K_5[, col_sel], max_mag_2K_5[, col_sel], col = color_2K)
points(max_doy_2K_6[, col_sel], max_mag_2K_6[, col_sel], col = color_2K)
points(max_doy_2K_7[, col_sel], max_mag_2K_7[, col_sel], col = color_2K)
points(max_doy_2K_8[, col_sel], max_mag_2K_8[, col_sel], col = color_2K)

points(max_doy_1p5K_1[, col_sel], max_mag_1p5K_1[, col_sel], col = color_1p5K)
points(max_doy_1p5K_2[, col_sel], max_mag_1p5K_2[, col_sel], col = color_1p5K)
points(max_doy_1p5K_3[, col_sel], max_mag_1p5K_3[, col_sel], col = color_1p5K)
points(max_doy_1p5K_4[, col_sel], max_mag_1p5K_4[, col_sel], col = color_1p5K)
points(max_doy_1p5K_5[, col_sel], max_mag_1p5K_5[, col_sel], col = color_1p5K)
points(max_doy_1p5K_6[, col_sel], max_mag_1p5K_6[, col_sel], col = color_1p5K)
points(max_doy_1p5K_7[, col_sel], max_mag_1p5K_7[, col_sel], col = color_1p5K)
points(max_doy_1p5K_8[, col_sel], max_mag_1p5K_8[, col_sel], col = color_1p5K)

points(max_doy_hist_1[, col_sel], max_mag_hist_1[, col_sel], col = color_hist)
points(max_doy_hist_2[, col_sel], max_mag_hist_2[, col_sel], col = color_hist)
points(max_doy_hist_3[, col_sel], max_mag_hist_3[, col_sel], col = color_hist)
points(max_doy_hist_4[, col_sel], max_mag_hist_4[, col_sel], col = color_hist)
points(max_doy_hist_5[, col_sel], max_mag_hist_5[, col_sel], col = color_hist)


