
### Processing Output

library(tidyverse)

inv_log <- function(x) exp(x)/(1 + exp(x))

# Read in Posterior Draws -------------------------------------------------

outs <- readRDS("outs.rds")

Pmod <- outs[[1]]
Tmod <- outs[[2]]

spp.table <- read.csv("spp.table.csv")

# make dataframe for global effects
global <- cbind.data.frame(med = c(median(Pmod$mu_prec), 
                                   median(Tmod$mu_temp)),
                           low = c(quantile(Pmod$mu_prec, 0.025), 
                                   quantile(Tmod$mu_temp, 0.025)),
                           up = c(quantile(Pmod$mu_prec, 0.975), 
                                  quantile(Tmod$mu_temp, 0.975)), 
                           var = c("prec", "temp"))

# make dataframe for species specific effects
species <- cbind.data.frame(species = unique(Pmod$s),
                            med = c(unname(tapply(Pmod$phi_prec, Pmod$s, median)), 
                                       unname(tapply(Tmod$phi_temp, Tmod$s, median))),
                            low = c(unname(tapply(Pmod$phi_prec, Pmod$s, quantile, 0.025)), 
                                    unname(tapply(Tmod$phi_temp, Tmod$s, quantile, 0.025))),
                            up = c(unname(tapply(Pmod$phi_prec, Pmod$s, quantile, 0.975)), 
                                    unname(tapply(Tmod$phi_temp, Tmod$s, quantile, 0.975))),
                            var = c(rep("prec", n_distinct(Pmod$s)), 
                                    rep("temp", n_distinct(Tmod$s))))

# replace species numbers with names
species$species <- spp.table$name[match(species$species, spp.table$number)]
# identify confidently non-zero effects
species$sig <- ifelse(species$low > 0 | species$up < 0, "yes", "no")


# Generate Table of Species Effects ---------------------------------------
effects_table <- cbind.data.frame(species = spp.table$name, 
                 median_temp_effect = as.numeric(tapply(Tmod$phi_temp, 
                                                        Tmod$s, median)), 
                 lower_temp_effect = as.numeric(tapply(Tmod$phi_temp, 
                                                        Tmod$s, quantile, 0.025)),
                 upper_temp_effect = as.numeric(tapply(Tmod$phi_temp, 
                                                       Tmod$s, quantile, 0.975)), 
                 median_prec_effect = as.numeric(tapply(Pmod$phi_prec, 
                                                        Pmod$s, median)), 
                 lower_prec_effect = as.numeric(tapply(Pmod$phi_prec, 
                                                       Pmod$s, quantile, 0.025)),
                 upper_prec_effect = as.numeric(tapply(Pmod$phi_prec, 
                                                       Pmod$s, quantile, 0.975)))


write.csv(effects_table, "species_effects.csv")


# Plotting Global Effects -------------------------------------------------

ggplot(data = global, aes(y = med, ymin = low, ymax = up, x = var)) + 
  geom_pointrange(fatten = 10, linewidth = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  xlab("Climate Variable") + ylab("Estimated Effect") + 
  theme_classic(base_size = 15) + 
  scale_x_discrete(labels = c("Precipitation", "Temperature"))

# Plotting Species Effects ------------------------------------------------

# specify order for y-axis
ordered.species <- species$species[order(species$med[which(species$var == "temp")])]
species$species <- factor(species$species, level = ordered.species)

# Precipitation Effects Graph
ggplot(data = species[which(species$var == "prec"),], 
       aes(y = species, x = med, xmin = low, xmax = up, 
           color = sig)) + 
  geom_pointrange() + 
  theme_classic(base_size = 15) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = .5), 
        legend.position = "none") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  ylab("Species") + xlab("Estimated Effect of Precipitation") + 
  scale_color_manual(values = c("grey60", "black"))

# Temperature Effects Graph
ggplot(data = species[which(species$var == "temp"),], 
       aes(y = species, x = med, xmin = low, xmax = up, 
           color = sig)) + 
  geom_pointrange() + 
  theme_classic(base_size = 15) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = .5), 
        legend.position = "none") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  ylab("Species") + xlab("Estimated Effect of Temperature") + 
  scale_color_manual(values = c("grey60", "black"))

# Get Informative Summaries -----------------------------------------------

# probability of directions
sum(Tmod$mu_temp < 0)/nrow(Tmod)
sum(species$sig[which(species$var=="temp")] == "yes")
min(tapply(Tmod$phi_temp, Tmod$s, function(x) sum(x<0)/length(x)))

sum(Pmod$mu_prec > 0)/nrow(Pmod)
sum(species$sig[which(species$var=="prec")] == "yes")
min(tapply(Pmod$phi_prec, Pmod$s, function(x) sum(x>0)/length(x)))

# bring in climate data to un-scale parameters & predictions

clim_dat <- read_table("climate R-Data.txt")
temps_scaled <- scale(clim_dat$Avg_temp[which(clim_dat$Year %in% 1985:2011)])
prec_scaled <- scale(clim_dat$Avg_rain[which(clim_dat$Year %in% 1985:2011)])

# What does a 1 degree increase in temp do?

delta_temp <- 1/attr(temps_scaled,"scaled:scale") # increase of 1 celsius on z-transformed scale

# % reduction in phi with 1 degree increase in temp
prob_red_temp <- ((inv_log(Tmod$mu_0) - 
                     inv_log(Tmod$mu_0 + Tmod$mu_temp*delta_temp))/
                    inv_log(Tmod$mu_0)) * 100

mean(prob_red_temp)

# What does a 10mm drop in prec do?

delta_prec <- -10/attr(prec_scaled,"scaled:scale") # decrease of 10 mm on z-transformed scale

# % reduction in phi with 10mm decrease in prec
prob_red_prec <- ((inv_log(Pmod$mu_0) - 
                     inv_log(Pmod$mu_0 + Pmod$mu_prec*delta_prec))/
                    inv_log(Pmod$mu_0)) * 100

median(prob_red_prec)
