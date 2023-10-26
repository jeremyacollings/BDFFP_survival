
# CHECKING MODEL OUTPUTS

library(tidyverse)
library(loo)

# bird_mod rds files were generated from the run_mods r scripts 
# they're sorta big, which is why I broke them up

### Bring in Data -----

mods1 <- readRDS('bird_mods1.rds')
mods2 <- readRDS('bird_mods2.rds')
mods3 <- readRDS('bird_mods3.rds')
mods4 <- readRDS('bird_mods4.rds')
mods5 <- readRDS('bird_mods5.rds')
mods6 <- readRDS('bird_mods6.rds')

# also need to bring in bird capture histories to get n.observed


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

n.observed <- c() # number of birds with an observation in the time period
for(i in sp){
  x <- bird_dat2[[i]]
  x[which(x == 3, arr.ind = TRUE)] <- 0
  n.observed <- c(n.observed,try(sum(rowSums(apply(x, 2, as.numeric), na.rm = TRUE) > 1)))
}

sp2 <- which(as.numeric(n.observed) >= 5)

### Checking Model Comparisons -----

# replace the first index with 1-5 to check different species
# the fourth item in the second dimension of the list is the 
## output of loo_compare()
# 1 = phi(.)p(.); 2 = phi(t)p(.); 3 = phi(.)p(t); 4 = phi(t)p(t); 
# 5 = phi(Prec)p(.); 6 = phi(Prec)p(t); 7 = phi(Temp)p(.); 8 = phi(Temp)p(t)

mods1[[1]][[4]] 

### Compiling Dataset for Viz -----

get_quant <- function(x, par, quant){
  unname(quantile(as.data.frame(x)[[par]], quant))
}

temp_mat0 <- matrix(NA, nrow = length(as.data.frame(mods1[[1]][[3]][[1]])[["phi_temp"]]), ncol = 29)
temp_mat <- matrix(NA, nrow = length(as.data.frame(mods1[[1]][[3]][[1]])[["phi_temp"]]), ncol = 29)
prec_mat0 <- matrix(NA, nrow = length(as.data.frame(mods1[[1]][[3]][[1]])[["phi_temp"]]), ncol = 29)
prec_mat <- matrix(NA, nrow = length(as.data.frame(mods1[[1]][[3]][[1]])[["phi_temp"]]), ncol = 29)

for(j in 1:5){
  temp_mat0[,j] <- as.data.frame(mods1[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,j] <- as.data.frame(mods1[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,j] <- as.data.frame(mods1[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,j] <- as.data.frame(mods1[[j]][[2]][[1]])[["phi_prec"]]
}
for(j in 1:5){
  temp_mat0[,5+j] <- as.data.frame(mods2[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,5+j] <- as.data.frame(mods2[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,5+j] <- as.data.frame(mods2[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,5+j] <- as.data.frame(mods2[[j]][[2]][[1]])[["phi_prec"]]
}
for(j in 1:5){
  temp_mat0[,10+j] <- as.data.frame(mods3[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,10+j] <- as.data.frame(mods3[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,10+j] <- as.data.frame(mods3[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,10+j] <- as.data.frame(mods3[[j]][[2]][[1]])[["phi_prec"]]
}
for(j in 1:5){
  temp_mat0[,15+j] <- as.data.frame(mods4[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,15+j] <- as.data.frame(mods4[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,15+j] <- as.data.frame(mods4[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,15+j] <- as.data.frame(mods4[[j]][[2]][[1]])[["phi_prec"]]
}
for(j in 1:5){
  temp_mat0[,20+j] <- as.data.frame(mods5[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,20+j] <- as.data.frame(mods5[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,20+j] <- as.data.frame(mods5[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,20+j] <- as.data.frame(mods5[[j]][[2]][[1]])[["phi_prec"]]
}
for(j in 1:4){
  temp_mat0[,25+j] <- as.data.frame(mods6[[j]][[3]][[1]])[["phi_0"]]
  temp_mat[,25+j] <- as.data.frame(mods6[[j]][[3]][[1]])[["phi_temp"]]
  prec_mat0[,25+j] <- as.data.frame(mods6[[j]][[2]][[1]])[["phi_0"]]
  prec_mat[,25+j] <- as.data.frame(mods6[[j]][[2]][[1]])[["phi_prec"]]
}

dat <- cbind.data.frame(species = c(names(mods1), names(mods2), 
                                    names(mods3), names(mods4), 
                                    names(mods5), names(mods6)), 
                        temp_med = apply(temp_mat, 2, median), 
                        temp_low = apply(temp_mat, 2, quantile, 0.025),
                        temp_up = apply(temp_mat, 2, quantile, 0.975), 
                        prec_med = apply(prec_mat, 2, median), 
                        prec_low = apply(prec_mat, 2, quantile, 0.025),
                        prec_up = apply(prec_mat, 2, quantile, 0.975), 
                        prec_sig = c(rep(TRUE, 12), FALSE, 
                                     rep(TRUE, 10), FALSE, rep(TRUE, 5)), 
                        temp_sig = c(rep(TRUE, 4), FALSE, rep(TRUE, 22), FALSE, TRUE), 
                        n = as.numeric(n.observed[sp2]))

dat$prec_med[which(dat$prec_sig == FALSE)] <- 0
dat$prec_low[which(dat$prec_sig == FALSE)] <- 0
dat$prec_up[which(dat$prec_sig == FALSE)] <- 0
dat$temp_med[which(dat$temp_sig == FALSE)] <- 0
dat$temp_low[which(dat$temp_sig == FALSE)] <- 0
dat$temp_up[which(dat$temp_sig == FALSE)] <- 0

dat$prec_sig2 <- dat$prec_low > 0
dat$temp_sig2 <- dat$temp_up < 0

### Forest Plots -----

# full temperature plot

ggplot(data = dat, aes(x = temp_med, y = reorder(species, temp_med), 
                       xmin = temp_low, xmax = temp_up, 
                       color = temp_sig2)) + 
  geom_point(size = 2) + geom_errorbar() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_classic(base_size = 15) + 
  xlab("Effect of Temperature on Survival") + 
  ylab("Species") + 
  scale_color_manual(name = element_blank(), values = c("#8E9CA4", "#253237")) + 
  theme(legend.position = "none", axis.text.y = element_text(size = 10))

# full precipitation plot

ggplot(data = dat, aes(x = prec_med, y = reorder(species, -prec_med), 
                       xmin = prec_low, xmax = prec_up, 
                       color = prec_sig2)) + 
  geom_point(size = 2) + geom_errorbar() + 
  geom_vline(xintercept = 0, linetype = "dashed") +  
  theme_classic(base_size = 15) + 
  xlab("Effect of Precipitation on Survival") + 
  ylab("Species") + 
  scale_color_manual(name = element_blank(), values = c("#8E9CA4", "#253237")) + 
  theme(legend.position = "none", axis.text.y = element_text(size = 10))

# truncated temperature plot

ggplot(data = dat[which(dat$n >= 8),], aes(x = temp_med, y = reorder(species, temp_med), 
                                           xmin = temp_low, xmax = temp_up, 
                                           color = temp_sig2)) + 
  geom_point(size = 2) + geom_errorbar() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_classic(base_size = 15) + 
  xlab("Effect of Temperature on Survival") + 
  ylab("Species") + 
  scale_color_manual(name = element_blank(), values = c("#8E9CA4", "#253237")) + 
  theme(legend.position = "none", axis.text.y = element_text(size = 10))

# truncated precipitation plot

ggplot(data = dat[which(dat$n >= 8),], aes(x = prec_med, y = reorder(species, -prec_med), 
                                           xmin = prec_low, xmax = prec_up, 
                                           color = prec_sig2)) + 
  geom_point(size = 2) + geom_errorbar() + 
  geom_vline(xintercept = 0, linetype = "dashed") +  
  theme_classic(base_size = 15) + 
  xlab("Effect of Precipitation on Survival") + 
  ylab("Species") + 
  scale_color_manual(name = element_blank(), values = c("#8E9CA4", "#253237")) + 
  theme(legend.position = "none", axis.text.y = element_text(size = 10))

