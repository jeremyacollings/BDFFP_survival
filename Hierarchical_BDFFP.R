
# Hierarchical CJS

library(rstan)
library(tidyverse)
library(tidybayes)
library(patchwork)
library(readxl)
library(loo)

options(mc.cores = parallel::detectCores())

# Bring in Data -----------------------------------------------------------

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

# bring in climate data
clim_dat <- read.table("climate R-Data.txt", header = TRUE)

# Prep Data for Stan ------------------------------------------------------

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

collapse.ch <- function(ch){
  ch.char = apply(ch, 1, function(x) paste(x, collapse = ","))
  
  ch.sum.out = t(sapply(strsplit(names(table(ch.char)), split = ","), as.numeric))
  fr.out = as.numeric(as.vector(table(ch.char)))
  
  return(list(ch.sum.out, fr.out))
}

get.first <- function(x) min(which(x != 0))

CHs1 <- CHs2 <- list()
sumfs <- sumFRs <- species <- c()
for(i in sp2){
  s = sp[i]
  CH0 <- as.matrix(apply(bird_dat2[[s]], 2, as.numeric))
  missing = which(CH0[1,] == 3)
  CH0[which(CH0 == 3, arr.ind = TRUE)] <- 0
  
  ch <- CH0
  ch2 <- ch[which(rowSums(ch) != 0), ]
  tmpCH = collapse.ch(ch2)[[1]]
  sumFR = collapse.ch(ch2)[[2]]
  
  sumCH = tmpCH
  sumCH2 = sumCH
  sumCH[sumCH[,] == 0] = 2
  
  CHs1[[i]] <- sumCH; CHs2[[i]] <- sumCH2
  
  sumfs <- c(sumfs,apply(tmpCH, 1, get.first))
  sumFRs <- c(sumFRs, sumFR)
  species <- c(species, rep(which(i == sp2), nrow(sumCH)))
}

CH_full1 <- do.call(rbind, CHs1)
CH_full2 <- do.call(rbind, CHs2)
NsumCH = nrow(CH_full1)         # number of capture histories 
n.occasions = ncol(CH_full1)    # number of sampling occasions

# get number of likelihoods

counter = 1;
for(j in 1:NsumCH){
  if (sumfs[j] >= n.occasions){
  }
  else{
    for(k in (sumfs[j]+1):(n.occasions)){
      for(l in 1:sumFRs[j]){
        counter = counter + 1;
      }
    }
  }
  
}

stan_dat <- list(NsumCH = NsumCH, 
                 n_occasions = n.occasions,
                 sumCH = CH_full1, 
                 CH = CH_full2,
                 sumf = sumfs, 
                 sumFR = sumFRs, 
                 n_missing = length(missing), 
                 missing = as.numeric(missing-1), 
                 n_observed = length(setdiff(1:ncol(CH_full1), missing))-1, 
                 observed = setdiff(1:ncol(CH_full1), missing)[-1]-1, 
                 temp = as.numeric(scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)])), 
                 prec = as.numeric(scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)])), 
                 n_lik = counter - 1, 
                 species = species, 
                 S = length(unique(species)))

# Fit Models --------------------------------------------------------------

phi.p. <- stan("Uncentered/phi(.)p(.)H.stan",
               data = stan_dat, 
               pars = c("log_lik", "p", "phi", "phi_0", "mu_0"))

phi.pt <- stan("Uncentered/phi(.)p(t)H.stan", 
               data = stan_dat,
               pars = c("log_lik", "p", "phi", "phi_0", "mu_0"))

phitp. <- stan("Uncentered/phi(t)p(.)H.stan", 
               data = stan_dat,
               pars = c("log_lik", "p", "phi", "phi_0",  "mu_0"))

phitpt <- stan("Uncentered/phi(t)p(t)H.stan", 
               data = stan_dat,
               pars = c("log_lik", "p", "phi", "phi_0", "mu_0"))

phiPRECp. <- stan("Uncentered/phi(prec)p(.)H.stan", 
                  data = stan_dat,
                  pars = c("log_lik", "p", "phi", "phi_0", "mu_0",
                           "mu_prec", "phi_prec"))

phiPRECpt <- stan("Uncentered/phi(prec)p(t)H.stan", 
                  data = stan_dat,
                  pars = c("log_lik", "p", "phi", "phi_0", "mu_0",
                           "mu_prec", "phi_prec"))

phiTEMPp. <- stan("Uncentered/phi(temp)p(.)H.stan", 
                  data = stan_dat,
                  pars = c("log_lik", "p", "phi", "phi_0", "mu_0",
                           "mu_temp", "phi_temp"))

phiTEMPpt <- stan("Uncentered/phi(temp)p(t)H.stan", 
                  data = stan_dat,
                  pars = c("log_lik", "p", "phi", "phi_0",  "mu_0",
                           "mu_temp", "phi_temp"))
  
# Model Comparison --------------------------------------------------------

loo_compare(loo(phi.p.), loo(phi.pt), 
            loo(phitp.), loo(phitpt), 
            loo(phiTEMPp.), loo(phiTEMPpt), 
            loo(phiPRECp.), loo(phiPRECpt))

# Export Model Output ----------------------------------------------------

phiPRECp.df <- spread_draws(model = phiPRECp., 
                         mu_0, p[s,t], phi_0[s], phi[s,t], 
                         mu_prec, phi_prec[s])

phiTEMPp.df <- spread_draws(model = phiTEMPp., 
                         mu_0, p[s,t], phi_0[s], phi[s,t],
                         mu_temp, phi_temp[s])

outs <- list(phiPRECp.df, phiTEMPp.df)

saveRDS(outs, "outs.rds")

# Calculation of ANODEV Statistic -----------------------------------------

prec_logs <- as.data.frame(phiPRECp.)[,which(grepl("log", names(as.data.frame(phiPRECp.))))]
temp_logs <- as.data.frame(phiTEMPp.)[,which(grepl("log", names(as.data.frame(phiTEMPp.))))]
null_logs <- as.data.frame(phi.p.)[,which(grepl("log", names(as.data.frame(phi.p.))))]
null_logs2 <- as.data.frame(phi.p.)[,which(grepl("log", names(as.data.frame(phitp.))))]

# calculate ANODEV stat from Program Mark book

perc_temp <- (-2*(rowSums(null_logs) + rowSums(temp_logs)))/
  (-2*(rowSums(null_logs) + rowSums(null_logs2)))

perc_prec <- (-2*(rowSums(null_logs) + rowSums(prec_logs)))/
  (-2*(rowSums(null_logs) + rowSums(null_logs2)))

quantile(perc_temp)
quantile(perc_prec)

write.csv(cbind.data.frame(temp = perc_temp, prec = perc_prec), 
          file = "ANODEV.csv")

