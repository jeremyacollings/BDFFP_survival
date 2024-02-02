library(loo)
library(rstanarm)

# Calculating some extra stuff

clim_dat <- read_table("climate R-Data.txt")
temps_scaled <- scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)])
prec_scaled <- scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)])

# function to import all sheets from excel file into a list
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- excel_sheets(filename)
  x <- lapply(sheets, function(X) read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

bird_dat <- read_excel_allsheets("Mark-recapture data BDFFP.xlsx")
sp <- names(bird_dat)[-length(bird_dat)]

# get just the years of interest
bird_dat2 <- list()
missing_years <- list() # just to check that the missing years is always the same
for(i in sp){
  bird_dat2[[i]] <- bird_dat[[i]][which(bird_dat[[i]][["rs"]] == 1501),]
  bird_dat2[[i]] <- bird_dat2[[i]][, which(names(bird_dat2[[i]]) %in% 1985:2012)]
  bird_dat2[[i]][which(bird_dat2[[i]] == ".", arr.ind = TRUE)] <- 3 # STAN won't like the dots
  missing_years[[i]] <- setdiff(1985:2012,as.numeric(names(bird_dat2[[i]]))) # which additional years are unsampled?
  for(j in missing_years[[i]]){
    bird_dat2[[i]][,as.character(j)] <- rep(3, nrow(bird_dat2[[i]]))
  }
  bird_dat2[[i]] <- bird_dat2[[i]][,order(names(bird_dat2[[i]]))] # get in the right order
}

# let's quickly retrieve the median estimates for everyone...

birds1 <- readRDS("bird_mods1.rds")
birds2 <- readRDS("bird_mods2.rds")
birds3 <- readRDS("bird_mods3.rds")
birds4 <- readRDS("bird_mods4.rds")
birds5 <- readRDS("bird_mods5.rds")
birds6 <- readRDS("bird_mods6.rds")
logs1 <- readRDS("log_liks1.rds")
logs2 <- readRDS("log_liks2.rds")
logs3 <- readRDS("log_liks3.rds")
logs4 <- readRDS("log_liks4.rds")
logs5 <- readRDS("log_liks5.rds")
logs6 <- readRDS("log_liks6.rds")

birds <- c(birds1, birds2, birds3, birds4, birds5, birds6)
logs <- c(logs1, logs2, logs3, logs4, logs5, logs6)

temp_points <- prec_points <- matrix(data = NA, nrow = 27, ncol = length(birds))
ests <- matrix(data = NA, nrow = length(birds), ncol = 13)
for(i in 1:length(birds)){
  temp_df <- as.data.frame(birds[[i]][["phiTEMPp."]][["model"]])
  prec_df <- as.data.frame(birds[[i]][["phiPRECp."]][["model"]])
  
  temp_points[,i] <- unname(apply(temp_df[,28:54], 2, median))
  prec_points[,i] <- unname(apply(prec_df[,28:54], 2, median))
  
  ests[i,1] <- names(birds)[[i]]
  
  ests[i,2] <- quantile(temp_df$phi_0, 0.025)
  ests[i,3] <- median(temp_df$phi_0)
  ests[i,4] <- quantile(temp_df$phi_0, 0.975)
  
  ests[i,5] <- quantile(temp_df$phi_temp, 0.025)
  ests[i,6] <- median(temp_df$phi_temp)
  ests[i,7] <- quantile(temp_df$phi_temp, 0.975)
  
  ests[i,8] <- quantile(prec_df$phi_0, 0.025)
  ests[i,9] <- median(prec_df$phi_0)
  ests[i,10] <- quantile(prec_df$phi_0, 0.975)
  
  ests[i,11] <- quantile(prec_df$phi_prec, 0.025)
  ests[i,12] <- median(prec_df$phi_prec)
  ests[i,13] <- quantile(prec_df$phi_prec, 0.975)
}

temp_sig <- prec_sig <- c()
for(i in 1:length(birds)){
  temp_birds1 <- which(rownames(birds[[i]][["full_comp"]]) == "model7") < 
    which(rownames(birds[[i]][["full_comp"]]) == "model1") &
    which(rownames(birds[[i]][["full_comp"]]) == "model7") <
    which(rownames(birds[[i]][["full_comp"]]) == "model2") &
    which(rownames(birds[[i]][["full_comp"]]) == "model7") <
    which(rownames(birds[[i]][["full_comp"]]) == "model3") &
    which(rownames(birds[[i]][["full_comp"]]) == "model7") <
    which(rownames(birds[[i]][["full_comp"]]) == "model4")
  temp_df <- as.data.frame(birds[[i]][["phiTEMPp."]][["model"]])
  temp_birds2 <- unname(quantile(temp_df$phi_temp, 0.025) > 0 | 
                    quantile(temp_df$phi_temp, 0.975) < 0)
  temp_sig <- c(temp_sig, temp_birds1&temp_birds2)
  
  prec_birds1 <- which(rownames(birds[[i]][["full_comp"]]) == "model5") < 
    which(rownames(birds[[i]][["full_comp"]]) == "model1") &
    which(rownames(birds[[i]][["full_comp"]]) == "model5") <
    which(rownames(birds[[i]][["full_comp"]]) == "model2") &
    which(rownames(birds[[i]][["full_comp"]]) == "model5") <
    which(rownames(birds[[i]][["full_comp"]]) == "model3") &
    which(rownames(birds[[i]][["full_comp"]]) == "model5") <
    which(rownames(birds[[i]][["full_comp"]]) == "model4")
  prec_df <- as.data.frame(birds[[i]][["phiPRECp."]][["model"]])
  prec_birds2 <- unname(quantile(prec_df$phi_prec, 0.025) > 0 | 
                         quantile(prec_df$phi_prec, 0.975) < 0)
  prec_sig <- c(prec_sig, prec_birds1&prec_birds2)
}

ests.df <- as.data.frame(ests)
names(ests.df) <- c("sp", "temp_int_low", "temp_int_med", "temp_int_up", 
                    "temp_slope_low", "temp_slope_med", "temp_slope_up", 
                    "prec_int_low", "prec_int_med", "prec_int_up", 
                    "prec_slope_low", "prec_slope_med", "prec_slope_up")

write.csv(ests.df, file = "estimates.csv")

# total number of captures...

n.catches <- c()
for(i in sp){
  x <- as.numeric(as.matrix(bird_dat2[[i]]))
  x[which(x == 3, arr.ind = TRUE)] <- 0
  n.catches <- c(n.catches, sum(x))
}
sum(n.catches)

# range of decreases with 1 degree increase in temp

delta_temp <- 1/attr(temps_scaled,"scaled:scale") # decrease of 1 celsius on z-transformed scale

# % reduction in phi with 1 degree increase in temp
prob_red_temp <- ((inv_log(as.numeric(ests.df$temp_int_med)) - 
                inv_log(as.numeric(ests.df$temp_int_med) + 
                          as.numeric(ests.df$temp_slope_med)*delta_temp))/
               inv_log(as.numeric(ests.df$temp_int_med))) * 100
range(prob_red_temp[which(temp_sig & n >= 7)])
# % reduction in survival odds with 1 degree increase in temp
odds_red_temp <- ((exp(as.numeric(ests.df$temp_int_med)) - 
                exp(as.numeric(ests.df$temp_int_med) + 
                      as.numeric(ests.df$temp_slope_med)*delta_temp))/
               exp(as.numeric(ests.df$temp_int_med))) * 100

odds_red_temp[which(temp_sig)]

# range of decreases with 10 mm dencrease in prec
delta_prec <- 10/attr(prec_scaled,"scaled:scale") # decrease of 1 celsius on z-transformed scale


# % reduction in phi with 10 mm dencrease in prec
prob_red_prec <- ((inv_log(as.numeric(ests.df$prec_int_med)) - 
                inv_log(as.numeric(ests.df$prec_int_med) + 
                          as.numeric(ests.df$prec_slope_med)*delta_prec))/
               inv_log(as.numeric(ests.df$prec_int_med))) * 100
range(prob_red_prec[which(prec_sig & n >= 7)])
# % reduction in survival odds with 10 mm dencrease in prec
odds_red <- ((exp(as.numeric(ests.df$prec_int_med)) - 
                exp(as.numeric(ests.df$prec_int_med) + 
                      as.numeric(ests.df$prec_slope_med)*delta_prec))/
               exp(as.numeric(ests.df$prec_int_med))) * 100

odds_red[which(prec_sig & n >= 7)]

# getting percent change in phi from coldest to warmest year

# most affected by temp
temp_sp <- which(ests.df$temp_slope_med == 
                   min(as.numeric(ests.df$temp_slope_med[which(n >=7)])))

points <- unname(apply(as.data.frame(birds[[temp_sp]][["phiTEMPp."]][["model"]])[,28:54], 
                2, median))
(points[which(temps_scaled == min(temps_scaled))] - 
  points[which(temps_scaled == max(temps_scaled))])/
  points[which(temps_scaled == min(temps_scaled))]

# most affected by prec
prec_sp <- which(ests.df$prec_slope_med == 
                   max(as.numeric(ests.df$prec_slope_med[which(n >=7)])))

points <- unname(apply(as.data.frame(birds[[prec_sp]][["phiPRECp."]][["model"]])[,28:54], 
                       2, median))
(points[which(prec_scaled == min(prec_scaled))] - 
    points[which(prec_scaled == max(prec_scaled))])/
  points[which(prec_scaled == min(prec_scaled))]


# likelihood calculations


perc_temps <- perc_prec <- matrix(data = NA, nrow = nrow(logs[[i]][["phi.p."]]), 
                                                         ncol = length(logs))
for(i in 1:length(logs)){
  perc_temps[,i] <- (-2*(rowSums(logs[[i]][["phi.p."]]) + 
                         rowSums(logs[[i]][["phiTEMPp."]])))/
                     (-2*(rowSums(logs[[i]][["phi.p."]]) + 
                            rowSums(logs[[i]][["phitp."]])))
  perc_prec[,i] <- (-2*(rowSums(logs[[i]][["phi.p."]]) + 
                         rowSums(logs[[i]][["phiPRECp."]])))/
                     (-2*(rowSums(logs[[i]][["phi.p."]]) + 
                            rowSums(logs[[i]][["phitp."]])))
}

range(apply(perc_temps, 2, median)[which(n >= 7)])
range(apply(perc_prec, 2, median)[which(n >= 7)])
