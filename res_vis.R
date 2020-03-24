###

#Analyze mhm model results

###

#set_up----

# devtools::install_github('ERottler/meltimr')
pacman::p_load(parallel, doParallel, zoo, zyp, alptempr, emdbook, scales)

run_dir <- "D:/nrc_user/rottler/mhm_run/6935053/"

bas_dir <- "U:/rhine_fut/R/"

#load functions
source(paste0(bas_dir, "mhm_ana/functs.R"))

sta_yea <- 1955
end_yea <- 2014

stopCluster(my_clust)

n_cores <- 25 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr))
registerDoParallel(my_clust)


#dis_ana----

dis_mhm <- read.table(paste0(run_dir, "output/daily_discharge.out"), header = T)
dis_mhm$date <- as.Date(strptime(paste0(dis_mhm$Day, "-", dis_mhm$Mon, "-", dis_mhm$Year), "%d-%m-%Y", tz="UTC"))

#Runoff seasonality

quants <- seq(0.01, 0.99, by = 0.01)

f_qvalu_obs <- function(quant_sel){dis_ana(disc = dis_mhm$Qobs_0006935053,
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

f_qvalu_sim <- function(quant_sel){dis_ana(disc = dis_mhm$Qsim_0006935053,
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


cols_max <- colorRampPalette(c(rep("grey98", 15), "lightgoldenrod2", "gold3", "goldenrod3", "orangered4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1, 1,2,3,4)], rep("lightcyan3", 1), rep("grey98", 15)))(100)
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
par(mar = c(2, 2, 0.5, 0.5))

vis_run(sta_day = "1955-01-01",
        end_day = "1957-12-31",
        do_legend = T)

vis_run(sta_day = "1982-01-01",
        end_day = "1984-12-31")

vis_run(sta_day = "2012-01-01",
        end_day = "2014-12-31")

dev.off()


#Plot: Discharge cumulative
pdf(paste0(bas_dir, "res_figs/runoff_cs.pdf"), width = 8, height = 4)

x_tics <- which(format(dis_mhm$date, "%m-%d") == "01-01")
x_labs <- 1951:2014

par(mar = c(3, 3, 0.5, 0.5))
plot(cumsum(dis_mhm$Qobs_0006935053), type = "n", ylab = "", xlab ="", axes = F)
lines(cumsum(dis_mhm$Qobs_0006935053), col = "black", lwd = 2)
lines(cumsum(dis_mhm$Qsim_0006935053), col = "red3", lwd = 2)
axis(2, mgp = c(3, 0.15, 0), tck = -0.02)
axis(1, at = x_tics, labels = x_labs)
mtext("Discharge cumulative [m³/s]", side = 2, line = 1.8)
legend("topleft", c("obs", "sim"), pch = 19, col = c("black", "red3"))
box()

dev.off()