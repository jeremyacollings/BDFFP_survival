
##### Surival Analyses for BDFFP Bird Data #####
# J.A. Collings 

library(tidyverse)
library(readxl)
library(rstan)
library(loo)
set.seed(6)

### Bringing in Data -----

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
clim_dat <- read_table("climate R-Data.txt")

### Processing Data -----

# get just the years of interest
bird_dat2 <- list()
missing_years <- list() # just to check that the missing years is always the same
for(i in sp){
  bird_dat2[[i]] <- bird_dat[[i]][, which(names(bird_dat[[i]]) %in% 1985:2012)]
  bird_dat2[[i]][which(bird_dat2[[i]] == ".", arr.ind = TRUE)] <- 0 # STAN won't like the dots
  missing_years[[i]] <- setdiff(1985:2012,as.numeric(names(bird_dat2[[i]]))) # which additional years are unsampled?
  for(j in missing_years[[i]]){
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

### Running Models -----

models <- c("CJS_...stan", "CJS_t..stan", "CJS_.t.stan", "CJS_tt.stan", 
            "CJS_P..stan", "CJS_W..stan", "CJS_P+W..stan", "CJS_PxW..stan", 
            "CJS_Pt.stan", "CJS_Wt.stan", "CJS_P+Wt.stan", "CJS_PxWt.stan")
# Note: For each time dependent p model, unsampled years have
# fixed p estimates at 0
rm.params <- c("chi", "mu") # I don't care about these parameters
sp <- names(birds_dat3)

# running peacemeal or I'll run out of memory
mods <- list()
loos <- list()
comps <- list()
for(s in sp[1:5]){
  stan_dat <- list(T = 28, I = nrow(birds_dat3[[s]]), y = as.matrix(apply(birds_dat3[[s]], 2, as.numeric)), 
                   temp = as.numeric(scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2012)])),
                   prec = as.numeric(scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2012)])), 
                   missing = which(1985:2012 %in% missing_years[[s]]), n_missing = length(missing_years[[s]]), 
                   observed = which(1985:2012 %in% setdiff(1985:2012, missing_years[[s]])), 
                   n_observed = length(setdiff(1985:2012, missing_years[[s]])))
  
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


