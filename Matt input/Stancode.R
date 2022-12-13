rm(list = ls())
#required packages
library(rstan)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Prject path
here::i_am("Stanmodel/Stancode.R")

#Data
growthdata <- read.delim("Data/Growthdata_pike.txt", 
                         header = T,stringsAsFactors = F,
                         sep = "\t", dec = ".")

#Data list for Stan
standata_growth = list(
  Radmm = growthdata$radmm, #numeric, otolith radius at age in mm (e.g., for age=2 radius = increment year 1 + increment year 2)
  Agei = growthdata$agei,   #integer, backcalculated age (age at radius measurement) in years
  N = nrow(growthdata),     #integer, number of observations
  n_fish = as.integer(max(growthdata$id3)), #integer, all fish numbers
  n_eco = as.integer(max(growthdata$ecotype)), #integer, number of ecotype groups (6)
  fish_id = as.integer(growthdata$id3)) #integer, individual id
fishtable = unique(growthdata[, c("id3", "ecotype")])
fishtable = fishtable[order(fishtable$id3), ]
standata_growth$eco_id_fish = fishtable$ecotype

#Stancode
stancode4 <- ("Data/Stancode/Growth_model_hierarchical4_ecotypes_3_normal.stan")

mod = stan_model(stancode4, model_name = "Individual pike growth")

#fitting
fit = sampling(mod, iter = 5000, data = standata_growth, 
           control = list(adapt_delta = 0.95, 
                          max_treedepth = 14))

#Investigate parameters of interest (Plausible range might be e.g., 
#LInf = 2.5 - 3.5, k = 0.1 - 0.5, tZero = -0.1 - -0.6)
print(fit,pars = c("LInf_eco","k_eco","tZero_eco","phi")) 
print(fit,pars = c("alpha","beta")) 
bayesplot::mcmc_hist(as.array(fit, pars = "phi"))
