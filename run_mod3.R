library(rstan)
library(tidyverse)
library(patchwork)
library(readxl)
library(loo)

##### Bringing in Data -----

# R2 function 
R2_func <- function(fit){
  fit.df <- as.data.frame(fit)
  ys <- as.matrix(fit.df[,which(grepl("y\\[", names(fit.df)))])
  y_hats <- as.matrix(fit.df[,which(grepl("y_hat", names(fit.df)))])
  Var_mu <- apply(y_hats, 1, var)
  Var_res <- rowSums(y_hats * (1-y_hats))/ncol(y_hats)
  Var_mu/(Var_mu+Var_res)
}

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

### Processing Data -----

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

##### Running STAN script -----

output <- list()
for(i in sp2[11:15]){
  s <- sp[i]
  
  CH <- as.matrix(apply(bird_dat2[[s]], 2, as.numeric))
  missing = which(CH[1,] == 3)
  CH[which(CH == 3, arr.ind = TRUE)] <- 0
  ch <- CH
  ch2 <- CH[which(rowSums(CH) != 0), ]
  collapse.ch <- function(ch){
    ch.char = apply(ch, 1, function(x) paste(x, collapse = ","))
    
    ch.sum.out = t(sapply(strsplit(names(table(ch.char)), split = ","), as.numeric))
    fr.out = as.numeric(as.vector(table(ch.char)))
    
    return(list(ch.sum.out, fr.out))
  }
  
  # format data for model fitting
  tmpCH = collapse.ch(ch2)[[1]]
  sumFR = collapse.ch(ch2)[[2]]
  
  # Create vector with occasion of marking
  get.first <- function(x) min(which(x != 0))
  sumf <- apply(tmpCH, 1, get.first)
  
  sumCH = tmpCH
  CH = sumCH
  sumCH[sumCH[,] == 0] = 2
  
  NsumCH = nrow(sumCH)         # number of capture histories 
  n.occasions = ncol(sumCH)    # number of sampling occasions
  
  # Catch (for the N versions)
  catch = colSums(CH)[2:18]
  
  # get number of likelihoods
  
  counter = 1;
  for(j in 1:NsumCH){
    if (sumf[j] >= n.occasions){
    }
    else{
      for(k in (sumf[j]+1):(n.occasions)){
        for(l in 1:sumFR[j]){
          counter = counter + 1;
        }
      }
    }
    
  }
  
  stan_dat <- list(NsumCH = NsumCH, 
                   n_occasions = n.occasions,
                   sumCH = sumCH, 
                   CH = CH,
                   sumf = sumf, 
                   sumFR = sumFR, 
                   n_missing = length(missing), 
                   missing = as.numeric(missing-1), 
                   n_observed = length(setdiff(1:ncol(CH), missing))-1, 
                   observed = setdiff(1:ncol(CH), missing)[-1]-1, 
                   temp = as.numeric(scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)])), 
                   prec = as.numeric(scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)])), 
                   n_lik = counter - 1)
  
  output[[s]][["phi.p."]][["model"]] <- phi.p. <- stan("stan_files/phi(.)p(.).stan",
                                                       data = stan_dat, 
                                                       cores = 4, chains = 4, 
                                                       pars = c("p", "phi", "log_lik", 
                                                                "y_hat"))
  output[[s]][["phi.p."]][["loo"]] <- phi.p._loo <- try(loo(phi.p.))
  phi.p._R2 <- R2_func(phi.p.)
  
  phitp. <- stan("stan_files/phi(t)p(.).stan", data = stan_dat, 
                 cores = 4, chains = 4, 
                 pars = c("p", "phi", "log_lik", 
                          "y_hat"))
  phitp._loo <- try(loo(phitp.))
  phitp._R2 <- R2_func(phitp.)
  
  phi.pt <- stan("stan_files/phi(.)p(t).stan", data = stan_dat, 
                 cores = 4, chains = 4, 
                 pars = c("p", "phi", "log_lik", 
                          "y_hat"))
  phi.pt_loo <- try(loo(phi.pt))
  phi.pt_R2 <- R2_func(phi.pt)
  
  phitpt <- stan("stan_files/phi(t)p(t).stan", data = stan_dat, 
                 cores = 4, chains = 4, 
                 pars = c("p", "phi", "log_lik", 
                          "y_hat"))
  phitpt_loo <- try(loo(phitpt))
  phitpt_R2 <- R2_func(phitpt)
  
  output[[s]][["phiPRECp."]][["model"]] <- phiPRECp. <- stan("stan_files/phi(prec)p(.).stan", data = stan_dat, 
                                                             cores = 4, chains = 4, 
                                                             pars = c("p", "phi", "phi_0", "phi_prec", "log_lik", 
                                                                      "y_hat"))
  output[[s]][["phiPRECp."]][["loo"]] <- phiPRECp._loo <- try(loo(phiPRECp.))
  phiPRECp._R2 <- R2_func(phiPRECp.)
  
  phiPRECpt <- stan("stan_files/phi(prec)p(t).stan", data = stan_dat, 
                    cores = 4, chains = 4, 
                    pars = c("p", "phi", "phi_0", "phi_prec", "log_lik", 
                             "y_hat"))
  phiPRECpt_loo <- try(loo(phiPRECpt))
  phiPRECpt_R2 <- R2_func(phiPRECpt)
  
  output[[s]][["phiTEMPp."]][["model"]] <- phiTEMPp. <- stan("stan_files/phi(temp)p(.).stan", data = stan_dat, 
                                                             cores = 4, chains = 4, 
                                                             pars = c("p", "phi", "phi_0", "phi_temp", "log_lik", 
                                                                      "y_hat"))
  output[[s]][["phiTEMPp."]][["loo"]] <- phiTEMPp._loo <- try(loo(phiTEMPp.))
  phiTEMPp._R2 <- R2_func(phiTEMPp.)
  
  phiTEMPpt <- stan("stan_files/phi(temp)p(t).stan", data = stan_dat, 
       cores = 4, chains = 4, 
       pars = c("p", "phi", "phi_0", "phi_temp", "log_lik", 
                "y_hat"))
  phiTEMPpt_loo <- try(loo(phiTEMPpt))
  phiTEMPpt_R2 <- R2_func(phiTEMPpt)
  
  output[[s]][["full_comp"]] <- try(loo_compare(phi.p._loo, phitp._loo, 
                                                phi.pt_loo, phitpt_loo, 
                                                phiPRECp._loo, phiPRECpt_loo, 
                                                phiTEMPp._loo, phiTEMPpt_loo))
  output[[s]][["R2s"]] <- list(phi.p._R2, phitp._R2, 
                            phi.pt_R2, phitpt_R2, 
                            phiPRECp._R2, phiPRECpt_R2, 
                            phiTEMPp._R2, phiTEMPpt_R2)
}

saveRDS(output, "bird_mods3.rds")
