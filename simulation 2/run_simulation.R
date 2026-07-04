##run simulation##

#load packages
library(GenKern)
library(MultiCure)
library(mstate)
library(mice)
library(flexsurv)
library(icmstate)

#load functions

source("generate_data.r")
source("run_scenario.r")
files <- list.files("../function", full.names = TRUE)
invisible(lapply(files, source))

#load scenarios for simulation 2
load("scenarios.RData")
full_factorial<-as.data.frame(full_factorial)

run_scenario(scen,n_sim=500,base_seed=12345,full_factorial) #scen can be chosen from 1 to 9.