### 

#Insert gauge Basel into mHM model set-up
#Erwin Rottler, Uni Potsdam, May 2020

###

run_dir <- "D:/nrc_user/rottler/mhm_run/6435060/"

#Projections
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

dem <- raster(paste0(run_dir, "input/morph/dem.asc"))
facc <- raster(paste0(run_dir, "input/morph/facc.asc"))
idgauges <- raster(paste0(run_dir, "input/morph/idgauges.asc"))

#Get meta data and measured time series for selected gauges (GRDC)

grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

coch_file <- paste0(grdc_dir, "6336050_Q_Day.Cmd.txt")
# base_file <- paste0(grdc_dir, "6935051_Q_Day.Cmd.txt")
gauge_ID <- 6336050

file_paths <- c(coch_file)

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

#gauge Cochem one row lower in facc
grdc_meta$latitude   <- grdc_meta$latitude - 0.05

gauges_84 <- SpatialPointsDataFrame(data.frame(lon = grdc_meta$longitude, 
                                               lat = grdc_meta$latitude), 
                                    data = data.frame(data =rep(-1, length(grdc_meta$longitude))), 
                                    proj4string =  crswgs84)
gauges    <- spTransform(gauges_84, epsg3035)

plot(facc, col = c("grey92", "blue2", "darkblue"), breaks = c(0, 5000, 1000000, 10000000), axes=F, legend = T, ylab = "", xlab = "", box = F)
points(gauges)

#Extract values from flow accumulation to verify that on river network
extract(facc, gauges)
extract(facc, gauges, buffer = 1000, df = T)

gauge_index <- which(values(facc) == extract(facc, gauges)); gauge_index
idgauges[gauge_index] <- gauge_ID
unique(values(idgauges))

#write new id gauges asc-file
writeRaster(idgauges, paste0(paste0(run_dir, "input/idgauges.asc")), overwrite = T)
#first value not -9999 but -9999.000000
#change manually in asc-file


#read discharge time series
# disc_base <- read_grdc(base_file)
disc_coch <- read_grdc(coch_file)

disc_yea <- format(disc_coch$date, "%Y")
disc_mon <- format(disc_coch$date, "%m")
disc_day <- format(disc_coch$date, "%d")
disc_hou <- rep("00", length(disc_coch$date))
disc_min <- rep("00", length(disc_coch$date))

disc_out <- cbind(disc_yea, disc_mon, disc_day, disc_hou, disc_min, disc_coch$value)

write.table(disc_out, paste0(run_dir, "input/", gauge_ID, ".day"), sep = "\t", row.names = F, col.names = F,
            quote = F)
#add header from existing gauging files manually

