###

#mHM model output analysis - functions

###

#trend statistics: basic function
dis_ana <- function(disc, date, start_year = 1950, end_year = 2010, method_analys, quant_in = .95,
                    method_quant = "gev", do_moving_average = T, window_width = 30, 
                    cover_thresh = 0.2, above_below_quant = "above", w_allign = "center", quant_annual = F,
                    break_day = 0, do_prob_beta_reg = F, beta_zero_one_thres = 0.2, rank_sel = 30,
                    weather_type = c(1, 2, 3, 4)){
  
  input_data_full <- data.frame(date = date, value = disc)
  
  #Clip selected time period
  input_data <- input_data_full[as.numeric(format(input_data_full$date,'%Y')) >= start_year, ]
  input_data <- input_data[as.numeric(format(input_data$date,'%Y')) <= end_year, ]
  
  #Fill possible gaps
  start_date <- as.POSIXct(strptime(paste0(start_year,"-01-01"), "%Y-%m-%d", tz="UTC"))
  end_date   <- as.POSIXct(strptime(paste0(end_year,"-12-31"),   "%Y-%m-%d", tz="UTC"))
  full_date  <- seq(start_date, end_date, by="day")
  
  input_data <- data.frame(dates  = full_date,
                           values = with(input_data, value[match(as.Date(full_date), as.Date(date))])
  )
  
  #Remove 29th of February
  input_data <- input_data[-which(format(input_data$date, "%m%d") == "0229"),]
  
  #Vector with the 365 days of the year
  days <- seq(as.Date('2014-01-01'), to=as.Date('2014-12-31'), by='days')
  days <- format(days,"%m-%d")
  
  #Order data by day
  data_day <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
  colnames(data_day) <- c("year", days)
  data_day[, 1] <- start_year:end_year
  
  if(break_day > 0){
    
    for(i in 0:(length(start_year:end_year) - 2)) {
      
      data_day[i+1, 2:366] <- input_data$values[(i*365 + 1 + break_day):((i+1)*365 + break_day)]
      
    }
    
    data_day <- data_day[-nrow(data_day),]
  }else{
    
    for(i in 0:(length(start_year:end_year)-1)) {
      
      data_day[i+1, 2:366] <- input_data$values[(i*365+1):((i+1)*365)]
      
    }
  }
  
  #Linear trend using Sen´s slope trend estimator
  f_sens_slope <- function(data_in, cover_thresh = 0.9){
    
    if(length(which(is.na(data_in))) / length(data_in) > (1-cover_thresh)){
      sens_slo <-  NA
    }else{
      time_step <- 1:length(data_in)
      sens_slo <- as.numeric(zyp.sen(data_in~time_step)$coefficients[2])
      #sens_slo <- as.numeric(zyp.trend.vector(data_in, method = "zhang", conf.intervals = F)[2])
    }
    return(sens_slo)
  }
  
  #Select quantile method
  
  if(method_quant == "empirical"){
    
    f_quan <- function(data_in){
      if((length(which(is.na(data_in))) / length(data_in)) > cover_thresh){#check cover threshold
        quant_out <-rep(NA, length(quant_in))}else{
          if(length(unique(data_in)) <= 1){
            quant_out <- rep(unique(data_in), length(quant_in))
          }else{
            quant_out <- quantile(data_in, probs = quant_in, type = 8, na.rm = T)}}
      
      return(as.numeric(quant_out))
      
    }
    
  }
  
  if(method_quant == "gev"){
    f_quan <- function(data_in){
      if((length(which(is.na(data_in))) / length(data_in)) > cover_thresh){#check cover threshold
        quant_out <- NA}else{
          if(length(unique(data_in)) <= 1){
            quant_out <- unique(data_in)
          }else{
            params <- lmr2par(data_in, type = "gev")
            if(is.null(params)){#check if parameters were calculated
              quant_out <- median(data_in)}else{
                quant_out <- par2qua(para = params, f = quant_in)}}}
      return(quant_out)
    }
  }
  
  if(method_quant == "gpd"){
    
    f_quan <- function(data_in){
      
      quant_out <- GPDquantile(data_in, probs = quant_in)
      
      return(quant_out)
      
    }
    
  }
  
  
  #Selecte analytical method
  
  if(method_analys == "mean"){
    f_mea_na_thres <- function(data_in){mea_na_thres(x = data_in, na_thres = 1 - cover_thresh)}
    res <- apply(data_day[,-1], 2, f_mea_na_thres)
  }
  
  if(method_analys == "median"){
    f_medi <- function(data_in){median(data_in, na.rm = T)}
    res <- apply(data_day[,-1], 2, f_medi)
  }
  
  if(method_analys == "sum"){
    f_sum <- function(data_in){sum(data_in, na.rm = T)}
    res <- apply(data_day[,-1], 2, f_sum)
  }
  
  if(method_analys == "sens_slope"){
    #Moving average filter
    input_data$ma <- rollapply(data = input_data$values, width = window_width,
                               FUN = mea_na_thres, align = "center", fill = NA)
    
    #Vector with the 365 days of the year
    days <- seq(as.Date('2014-01-01'), to=as.Date('2014-12-31'), by='days')
    days <- format(days,"%m-%d")
    
    #Order data by day
    data_day <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
    colnames(data_day) <- c("year", days)
    data_day[ ,1] <- start_year:end_year
    
    for(i in 0:(length(start_year:end_year)-1)) {
      
      data_day[i+1, 2:366] <- input_data$ma[(i*365+1):((i+1)*365)]
      
    }
    
    #Trends only calculated when at least 50 % of input not NA or 0
    res <- apply(data_day[,-1], 2, f_sens_slope)
    
    
  }
  
  if(method_analys == "quantile"){
    
    if(quant_annual){
      res <- apply(data_day[,-1], 1, f_quan)
    }else{
      res <- apply(data_day[,-1], 2, f_quan)
    }
  }
  
  if(method_analys == "mov_quant_trend"){
    
    
    ns <- 1:(length(input_data$values) - 29) #interger to moving window over time series
    
    mov_quants <- function(int){
      
      if(length(which(is.na(input_data$values[int:(int+window_width-1)]))) > (window_width/2)){
        
        mov_quant_out <- rep(NA, length(quant_in))
        
      }else{
        mov_quant_out <- f_quan(input_data$values[int:(int+window_width-1)])
        
      }
      
      return(mov_quant_out)
      
    }
    
    mov_quants_vals <- pbsapply(ns, mov_quants, simplify = T)
    
    #add NA columns beginning and end to 'center' moving window
    
    na_cols_sta <- matrix(NA, nrow = nrow(mov_quants_vals), ncol = 14)
    na_cols_end <- matrix(NA, nrow = nrow(mov_quants_vals), ncol = 15)
    mov_quants_vals <- cbind(na_cols_sta, mov_quants_vals, na_cols_end)
    
    #Calculate trends of quantiles
    
    f_quant_wind_slo <- function(row_sel){
      
      data_quant_sel <- mov_quants_vals[row_sel, ]
      
      #Order data by day
      data_day <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
      colnames(data_day) <- c("year", days)
      data_day[ ,1] <- start_year:end_year
      
      for(i in 0:(length(start_year:end_year)-1)) {
        
        data_day[i+1, 2:366] <- data_quant_sel[(i*365+1):((i+1)*365)]
        
      }
      
      res_quant <- apply(data_day[, -1], 2, f_sens_slope) * 10 # per decade
      
      return(res_quant)
      
    }
    
    quants_int <- 1:nrow(mov_quants_vals)
    
    res <- pbsapply(quants_int, f_quant_wind_slo, simplify = T)
    
    
    
    # #Moving quantiles
    # input_data$mq <- rollapply(data = input_data$values, width = window_width,
    #                            FUN = f_quan, align = "center", fill = NA)
    # 
    # #Order data by day
    # data_day <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
    # colnames(data_day) <- c("year", days)
    # data_day[ ,1] <- start_year:end_year
    # 
    # for(i in 0:(length(start_year:end_year)-1)) {
    #   
    #   data_day[i+1, 2:366] <- input_data$mq[(i*365+1):((i+1)*365)]
    #   
    # }
    # 
    # #Calculate trends of quantiles
    # f_sens_slope <- function(data_in, cover_thresh = 0.2){
    #   if(length(which(is.na(data_in))) / length(data_in) > (1-cover_thresh)){
    #     sens_slo <-  NA
    #   }else{
    #     sens_slo <- as.numeric(zyp.trend.vector(data_in, method = "zhang", conf.intervals = F)[2])
    #   }
    #   return(sens_slo)
    # }
    # res <- apply(data_day[, -1], 2, f_sens_slope)
    
  }
  
  return(res)
  
}

#Plot: Discharge image
dis_image <- function(data_plot, cols, breaks, header, lab_unit, do_cont = T){
  
  par(mar = c(1.8, 1.8, 1.8, 0.2))
  par(family = "serif")
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  ytiks      <- seq(10, 90, by =  10)
  ylabs      <- seq(10, 90, by =  10)
  
  image(x = 1:365,
        y = 1:ncol(data_plot),
        z = data_plot, col = cols, breaks = breaks,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
  axis(2, at = ytiks, labels = ylabs/100, mgp=c(3, 0.15, 0), tck = -0.02, cex.axis = 1.2)
  mtext(header, side = 3, line = 0.3, cex = 1.1, adj = 0.0)
  mtext(lab_unit, side = 3, line = 0.2, cex = 1.0, adj = 1.0)
  box()
  
  if(do_cont){
    
    probs_iso <- c(0.1, 0.5, 0.9)
    iso_def <- quantile(data_plot, probs = probs_iso, type = 8, na.rm = T)
    
    contour(x = 1:365,
            y = 1:ncol(data_plot),
            z = as.matrix(data_plot),
            levels = round(iso_def, 0),
            add = T,
            lwd = 0.7,
            labcex = 0.7)
    
  }
  
  
  par(mar = c(1.8, 0.2, 1.8, 3.0))
  
  alptempr::image_scale(as.matrix(data_plot), col = cols, breaks = breaks, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
  box()
  
  
}

#Plot: Discharge time series
vis_run <- function(sta_day, end_day, do_legend = F){
  
  ind_sta <- which(dis_mhm$date == sta_day)
  ind_end <- which(dis_mhm$date == end_day)
  
  col_obs <- alpha("black", alpha = 0.6)
  col_sim <- alpha("red3", alpha = 0.6)
  
  ylims <- range(c(dis_mhm$Qobs_0006435060[ind_sta:ind_end], dis_mhm$Qsim_0006435060[ind_sta:ind_end]))
  plot(dis_mhm$Qobs_0006435060[ind_sta:ind_end], type = "n", ylim = ylims, axes = F)
  lines(dis_mhm$Qobs_0006435060[ind_sta:ind_end], col = col_obs, lwd = 1.5)
  lines(dis_mhm$Qsim_0006435060[ind_sta:ind_end], col = col_sim, lwd = 1.5)
  axis(1, at = c(1, length(ind_sta:ind_end)), labels = c(sta_day, end_day), mgp = c(3, 0.15, 0), tck = -0.02)
  axis(2, mgp = c(3, 0.15, 0), tck = -0.02)
  mtext("Discharge [m³/s]", side = 2, line = 1.7, cex = 0.7)
  if(do_legend){
    legend("topleft", c("obs", "sim"), pch = 19, col = c("black", "red3"))
  }
  box()
  
}

#Values to colors for image plot
val2col <- function(val_in, dat_ref, do_log = F, do_bicol = T, col_na = "white", virid_dir = -1){
  
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
    if(virid_dir == -1){
      my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
    }
    
    if(virid_dir == 1){
      my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = 1)))(200))
    }
    
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

#Get discharge from nc-files
