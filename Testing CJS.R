# gonna try to work with M. longipennis because it 
# is supposedly the most sensitive to rain and temp

library(RMark)
library(rstan)
library(tidyverse)
library(loo)
library(readxl)

# bring in all data

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- excel_sheets(filename)
  x <- lapply(sheets, function(X) read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

bird_dat <- read_excel_allsheets("Mark-recapture data BDFFP.xlsx")
sp <- names(bird_dat)[-length(bird_dat)]

# just the years of interest... ie 1985:2012

bird_dat2 <- list()
for(i in sp){
  bird_dat2[[i]] <- bird_dat[[i]][, which(names(bird_dat[[i]]) %in% 1985:2012)]
  bird_dat2[[i]][which(bird_dat2[[i]] == ".", arr.ind = TRUE)] <- 0 # FIGURE OUT UNSAMPLED
  missing_years <- setdiff(1985:2012,as.numeric(names(bird_dat2[[i]]))) # which additional years are unsampled?
  for(j in missing_years){
    bird_dat2[[i]][,as.character(j)] <- rep(0, nrow(bird_dat2[[i]]))
  }
  bird_dat2[[i]] <- bird_dat2[[i]][,order(names(bird_dat2[[i]]))] # get in the right order
}

n.observed <- c() # number of birds with an observation in the time period
for(i in sp){
  x <- bird_dat2[[i]]
  n.observed <- c(n.observed,sum(rowSums(apply(x, 2, as.numeric), na.rm = TRUE) > 1))
}

birds_dat3 <- bird_dat2[which(n.observed > 25)] # get rid of sp with less than 25 birds observed in time period

# bring in climate data
clim_dat <- read_table("climate R-Data.txt")

models <- c("CJS_...stan", "CJS_t..stan", "CJS_.t.stan", "CJS_tt.stan", 
            "CJS_P..stan", "CJS_W..stan", "CJS_P+W..stan", "CJS_PxW..stan", 
            "CJS_Pt.stan", "CJS_Wt.stan", "CJS_P+Wt.stan", "CJS_PxWt.stan")
#models <- c("CJS_...stan", "CJS_t..stan")
rm.params <- c("chi", "mu")

sp <- names(birds_dat3)
#sp <- names(birds_dat3)[1:2]
mods <- list()
loos <- list()
comps <- list()
for(s in sp[1:5]){
  stan_dat <- list(T = 28, I = nrow(birds_dat3[[s]]), y = as.matrix(apply(birds_dat3[[s]], 2, as.numeric)), 
                   temp = as.numeric(scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2012)])),
                   prec = as.numeric(scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2012)])))
  
  for(i in models){
    mods[[s]][[i]] <- stan(i,
                             data = stan_dat,
                             pars = rm.params,
                             include = FALSE,
                             control = list(adapt_delta = .85),
                             chains = 4, iter = 2000) 
    loos[[s]][[i]] <- loo(mods[[s]][[i]])
  }
  
  comps[[s]] <- loo_compare(loos[[s]])
}

saveRDS(mods, "mods11-15.RDS")
saveRDS(loos, "loos11-15.RDS")
saveRDS(comps, "comps11-15.RDS")


mods <- list()
loos <- list()
comps <- list()
for(s in sp[16:20]){
  stan_dat <- list(T = 28, I = nrow(birds_dat3[[s]]), y = as.matrix(apply(birds_dat3[[s]], 2, as.numeric)), 
                   temp = as.numeric(scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2012)])),
                   prec = as.numeric(scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2012)])))
  
  for(i in models){
    mods[[s]][[i]] <- stan(i,
                           data = stan_dat,
                           pars = rm.params,
                           include = FALSE,
                           control = list(adapt_delta = .85),
                           chains = 4, iter = 2000) 
    loos[[s]][[i]] <- loo(mods[[s]][[i]])
  }
  
  comps[[s]] <- loo_compare(loos[[s]])
}
# This crashed my computer... let's break it up... maybe five at a time???

# bind them together

mods_full <- c(read_rds("mods1-5.RDS"), read_rds("mods6-10.RDS"))
loos_full <- c(read_rds("loos1-5.RDS"), read_rds("loos6-10.RDS"))
comps_full <- c(read_rds("comps1-5.RDS"), read_rds("comps6-10.RDS"))

for(i in 1:10){
  print(rownames(comps_full[[i]])[1:2]) # print best two model names
}

# okay so... the best model is always a time variant survival model
# second best is often the prec and warming interaction model!
# if not thought, the second best is the null
# but how much better??

for(i in 1:10){
  print(comps_full[[i]][2,1:2]) # how much better is the best model than the second best
}

# only one of these distinctions is negligable, and its one of the null models

comps_full[[4]] # there's pretty similar model performance for most models 
# with time invariant detection probability

expit <- function(x) 1/(1+exp(-x))
unscale <- function(x, mu, sd) mu+x*sd
plots <- list()
for(i in 1:10){
  output <- as.data.frame(mods_full[[i]][["CJS_P+W..stan"]])
  sp <- names(mods_full)[i]
  plot_dat <- cbind.data.frame(surv = apply(output[,1:27], 2, function(x) median(expit(x))), 
                               temp = clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)], 
                               rain = clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)], 
                               #inter = clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)]*
                               #  clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)], 
                               temp_scaled = scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)]), 
                               rain_scaled = scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)]))
                               #inter_scaled = scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)])*
                               #  scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)]))
  
  plots[[sp]][['temp']] <- ggplot(data = plot_dat, aes(x = temp, y = surv)) + geom_point(size = 2) + 
    theme_classic(base_size = 15) + xlab("Temperature (C)") + 
    ylab("Estimated Survival Probability") + ggtitle(sp)
  
  plots[[sp]][['rain']] <- ggplot(data = plot_dat, aes(x = rain, y = surv)) + geom_point(size = 2) + 
    theme_classic(base_size = 15) + xlab("Precipitation") + 
    ylab("Estimated Survival Probability") + ggtitle(sp)
  
  # plots[[sp]][['inter']] <- ggplot(data = plot_dat, aes(x = inter, y = surv)) + geom_point(size = 2) + 
  #   theme_classic(base_size = 15) + xlab("Temp * Prec") + 
  #   ylab("Estimated Survival Probability") + ggtitle(sp)
  # 
  plot_dat2 <- cbind.data.frame(par = c("Temp", "Rain"), #, "Interaction"), 
                                med = c(median(output$phi_temp), 
                                        median(output$phi_prec)),#, 
                                        #median(output$phi_inter)),
                                low = c(quantile(output$phi_temp, 0.025), 
                                        quantile(output$phi_prec, 0.025)),#, 
                                        #quantile(output$phi_inter, 0.025)), 
                                up = c(quantile(output$phi_temp, 0.975), 
                                       quantile(output$phi_prec, 0.975))) #, 
                                       #quantile(output$phi_inter, 0.975)))
  
  plots[[sp]][["est"]] <- ggplot(data = plot_dat2, aes(x = par, y = med)) + 
    geom_point(size = 3) + 
    geom_errorbar(aes(ymin = low, ymax = up), width = .25) + 
    theme_classic(base_size = 15) + 
    xlab("Predictor") + ylab("Estimate") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle(sp)
}

for(i in names(plots)){
  plot(plots[[i]][["est"]])
}
                              