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
  fish_id = as.integer(growthdata$id3), #integer, individual id
  eco_id = as.integer(growthdata$ecotype)) #integer, ecotype id 1 - 6

#Stancode
stancode4 <- ("Data/Stancode/Growth_model_hierarchical4_ecotypes_3_normal.stan")

#fitting
fit = stan(stancode4, iter = 1000, data = standata_growth, 
           model_name = "Individual pike growth",
           control = list(adapt_delta = 0.9, 
                          max_treedepth = 12))

#Investigate parameters of interest (Plausible range might be e.g., 
#LInf = 2.5 - 3.5, k = 0.1 - 0.5, tZero = -0.1 - -0.6)
print(fit,pars = c("LInf_eco","k_eco","tZero_eco","phi")) 

