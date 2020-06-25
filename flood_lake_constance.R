###

#Influence of Lake Constance on flood wave propagation
#Erwin Rottler, University of Potsdam

###

bas_dir <- "U:/rhine_fut/R/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"
data_dir <- "D:/nrc_user/rottler/"

disc_diep <- meltimr::read_grdc(paste0(grdc_dir, "6935500_Q_Day.Cmd.txt"))
disc_neuh <- meltimr::read_grdc(paste0(grdc_dir, "6935055_Q_Day.Cmd.txt"))
disc_reki <- meltimr::read_grdc(paste0(grdc_dir, "6935054_Q_Day.Cmd.txt"))
disc_base <- meltimr::read_grdc(paste0(grdc_dir, "6935051_Q_Day.Cmd.txt"))

diep_day <- meltimr::ord_day(data_in = disc_diep$value,
                             date = disc_diep$date,
                             start_y = 1919,
                             end_y = 2016,
                             break_day = 0)

neuh_day <- meltimr::ord_day(data_in = disc_neuh$value,
                             date = disc_neuh$date,
                             start_y = 1919,
                             end_y = 2016,
                             break_day = 0)

reki_day <- meltimr::ord_day(data_in = disc_reki$value,
                             date = disc_reki$date,
                             start_y = 1919,
                             end_y = 2016,
                             break_day = 0)

base_day <- meltimr::ord_day(data_in = disc_base$value,
                             date = disc_base$date,
                             start_y = 1919,
                             end_y = 2016,
                             break_day = 0)

pdf(paste0(bas_dir,"res_figs/infl_const.pdf"), width = 8, height = 11)

comp_disc_plot <- function(year_sel){
  
  par(mar = c(2, 3.5, 2, 1))
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  year_all <- 1919:2016
  year_ind <- which(year_all == year_sel)
  
  ylims <- c(0, meltimr::max_na(c(diep_day[year_ind, ], neuh_day[year_ind, ], reki_day[year_ind, ], base_day[year_ind, ])))
  
  col_base <- viridis::viridis(100)[22]
  col_reki <- viridis::viridis(100)[45]
  col_neuh <- "grey25"
  col_diep <- "gold3"
  
  plot(x = c(1:365), diep_day[year_ind, ], type = "n", ylim = ylims, ylab = "", xlab = "", 
       axes = F, xaxs = "i", yaxs = "i")
  abline(v = x_axis_tic, lty = "dotted", col = "grey55", lwd = 0.8)
  polygon(x = c(1:365, 365:1), y = c(base_day[year_ind, ], rev(reki_day[year_ind, ])),
          col = scales::alpha(col_base, alpha = 0.7), border = F)
  polygon(x = c(1:365, 365:1), y = c(reki_day[year_ind, ], rev(neuh_day[year_ind, ])),
          col = scales::alpha(col_reki, alpha = 0.7), border = F)
  polygon(x = c(1:365, 365:1), y = c(neuh_day[year_ind, ], rev(diep_day[year_ind, ])),
          col = scales::alpha(col_neuh, alpha = 0.7), border = F)
  polygon(x = c(1:365, 365:1) , y = c(diep_day[year_ind, ], rep(0, 365)),
          col = scales::alpha(col_diep, alpha = 0.7), border = F)
  
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.04)#plot ticks
  axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S","O","N","D"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.08, 0), cex.axis = 1.3)#plot labels
  axis(2, mgp=c(3, 0.08, 0), tck = -0.01, cex.axis = 1.3)
  box()
  mtext(expression(paste("Discharge", " [m"^"3", "s"^"-1", "]")), side = 2, line = 1.4, cex = 1.1, adj = 0.5)
  mtext(year_sel, side = 3, line = 0.04, adj = 0.5, cex = 1.1)
  legend("topleft", c("Basel", "Rekingen", "Neuhausen", "Diepoldsau"), col = c(col_base, col_reki, col_neuh, col_diep),
         pch = 19, cex = 1.1, bg = "white")
  
}

par(mfrow = c(4, 1))
      
comp_disc_plot(1930)
comp_disc_plot(1940)
comp_disc_plot(1990)
comp_disc_plot(2000)

dev.off()


#map_over----

dem = raster::raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_500.tif"))

rivers <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Large_rivers.shp"))
tribus <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Other_large_rivers_and_tributaries.shp"))
lakes  <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Large_lakes.shp"))

rhin_riv <- rivers[rivers@data$NAME == "Rhine",]
bodensee <- lakes[which(lakes@data$WPLKNM == "BODENSEE" ),]

basins <-  rgdal::readOGR(dsn = paste0(data_dir, "basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp"), encoding = "UTF8")

basin_base <- sp::spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = raster::crs(dem, asText = T))
basin_neuh <- sp::spTransform(basins[basins@data$Ort == "Neuhausen, Flurlingerbrücke",], CRS = raster::crs(dem, asText = T))
basin_reki <- sp::spTransform(basins[basins@data$Ort == "Rekingen",], CRS = raster::crs(dem, asText = T))
basin_diep <- sp::spTransform(basins[basins@data$Ort == "Diepoldsau, Rietbrücke",], CRS = raster::crs(dem, asText = T))

#7.6167, 47.5594 Basel
#8.33, 47.5704 Rekingen
#8.6259, 47.6823 Neuhausen
#9.6409, 47.3831 Diepoldsau
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gauges_84 <-  sp::SpatialPoints(data.frame(lon = c(7.6167, 8.33, 8.6259, 9.6409),
                                           lat = c(47.5594, 47.5704, 47.6823, 47.3831)),
                                proj4string =  crswgs84)
gauges    <- sp::spTransform(gauges_84, CRS = raster::crs(basin_base, asText = T))

#corp DEM sub-basin area
my_ext <- raster::extent(basin_base)
my_ext_buf <- my_ext + c(-20000, +30000, -10000, +13000) #xmin, xmax, ymin, ymax

my_box <- as(my_ext_buf, 'SpatialPolygons')
dem_cro_swiss <- raster::crop(dem, raster::extent(my_box))
dem_sub_swiss <- raster::mask(dem_cro_swiss, my_box)

pdf(paste0(bas_dir,"res_figs/map_over.pdf"), width = 10, height = 6)

par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
par(family = "serif")

col_base <- viridis::viridis(100)[22]
col_reki <- viridis::viridis(100)[45]
col_neuh <- "grey25"
col_diep <- "gold3"

raster::image(dem_sub_swiss, axes = F, legend = F,  col = grDevices::colorRampPalette(c("white", "black"))(200), box = F)
raster::plot(basin_base, add =T, lwd = 0.025, col = scales::alpha(col_base, alpha = 0.45))
raster::plot(basin_reki, add =T, lwd = 0.025, col = scales::alpha(col_reki, alpha = 0.45))
raster::plot(basin_neuh, add =T, lwd = 0.025, col = scales::alpha(col_neuh, alpha = 0.65))
raster::plot(basin_diep, add =T, lwd = 0.025, col = scales::alpha(col_diep, alpha = 0.55))

raster::plot(rhin_riv, add = T, col = "blue4")
# plot(aare_riv, add = T, col = "blue4")
# plot(brienzersee, add = T, col = "blue4")
# plot(thunersee, add = T, col = "blue4")
# plot(bielersee, add = T, col = "blue4")
raster::plot(bodensee, add = T, col = "blue4")
raster::plot(gauges, add = T, pch = 23, cex = 1.7,
             bg = scales::alpha("grey55", alpha = 1.0))
raster::plot(gauges, add = T, pch = 19, cex = 0.5)
lab_mov <- 7000
lab_pos_1 <- c(0, 0, 0, 0)
lab_pos_2 <- c(-lab_mov, -lab_mov, lab_mov, -lab_mov)
text(gauges@coords[, 1]+lab_pos_1, gauges@coords[, 2]+lab_pos_2, 
     labels = c("Basel", "Rekingen", "Neuhausen", "Diepoldsau"), col = "black", cex = 1.2)
prettymapr::addscalebar(plotunit = "m", widthhint = 0.2, htin = 0.15, pos = "topleft",
                        padin = c(0.2, 0.2))

dev.off()








f_doy_max <- function(data_in){
  
  doy_max <- which(data_in == meltimr::max_na(data_in))[1]
  
  return(doy_max)
  
}


diep_doy <- apply(diep_day, 1, f_doy_max)
neuh_doy <- apply(neuh_day, 1, f_doy_max)

plot(diep_doy - neuh_doy)
abline(h = 0)
hist(diep_doy-neuh_doy, nclass =20)

