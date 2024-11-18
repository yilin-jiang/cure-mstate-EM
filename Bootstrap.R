library(mstate)
library(survival)
library(stringr)
library(doParallel)
library(foreach)
source("create_datalonglong.R")
source("update_pi.R")
source("calculate_hazard.R")
source("add_likelihood.R")
source("calculate_weight_mstate.R")
source("create_datalogit.R")
source("calculate_loglik_ob.R")
source("bootstrap_procedure.R")

######## SETTING #######
data("ebmt4")
data=ebmt4
tmat=transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                       c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
time=c("rec", "ae","recae", "rel", "srv")
status=c("rec.s", "ae.s", "recae.s","rel.s", "srv.s")
var=c("match", "proph", "year", "agecl")
subjectid='id'
tran_cure=c("12","13","16","24","26","34","36","46")
beta_fit=readRDS("cox_mstate.rda")
alpha_fit=readRDS("fit_alpha.rda")
data_coxfit=read.table("data_longlong.csv",sep=",",head=TRUE)

###########################
#### PARALLEL COMPUTE #####
N.boot=1000

# Number of cores to use for parallel processing
n_cores = detectCores()

# Create a parallel backend using doParallel
cl = makeCluster(n_cores)
registerDoParallel(cl)

bootstrap_results=foreach(boot_id = 1:N.boot,.export = c(".GlobalEnv"), .combine="rbind",
                          .packages=c("survival","stringr","mstate")) %dopar% {
                boot_pro(boot_id = boot_id,
                         data=ebmt4,
                         data_coxfit=data_coxfit,
                         time=c("rec", "ae","recae", "rel", "srv"),
                         status=c("rec.s", "ae.s", "recae.s","rel.s", "srv.s"),
                         tmat=tmat,
                         var=c("match", "proph", "year", "agecl"),
                         tran_cure=c("12","13","16","24","26","34","36","46"),#from-to pairs possible for cured subjects
                         alpha_fit=alpha_fit,
                         beta_fit=beta_fit,
                         subjectid="id",
                         max.iter=80)
                          }

# Stop the parallel backend
stopCluster(cl)