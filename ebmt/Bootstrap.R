# Load required packages
library(mstate)
library(survival)
# Load required functions
files <- list.files("../function", full.names = TRUE)
invisible(lapply(files, source))
# Load results from the core run
result_elf<-readRDS("ebmt results/result_elf.rda")

# SETTING 
data("ebmt4")
data<-ebmt4
tmat<-transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                        c(6), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
N.boot1<-1
N.boot2<-1000
  
master_seed<-2026

results_boot <- vector("list", length = N.boot2 - N.boot1 + 1)
names(results_boot) <- paste0("iter_", N.boot1:N.boot2)
failed_runs <- data.frame(iteration = integer(), error = character(), 
                          stringsAsFactors = FALSE)

# Bootstrap procedure
for (i in N.boot1:N.boot2){
  data_boot<-gendataboot(data=data,seed=master_seed+i)
  data_long_boot<- msprep(data = data_boot, trans = tmat, 
                        time = c(NA, c("rec", "ae","recae", "rel", "srv")), 
                         status = c(NA,c("rec.s", "ae.s", "recae.s","rel.s", "srv.s")), 
                         keep = c("match", "proph", "year", "agecl"))
 
 boot_result=tryCatch(
    {
      EM_mstate_cure( data_long0 = data_long_boot,
      id           = "id",
      covar_cure   = c("year","agecl","proph","match"),
      covar_tran   = rep(list(c("year","agecl","proph","match")), 13),
      tran_state   = c("12","13","15","16","24","25","26","34","35","36","45","46","56"),
      tran_noncure = c("1","2","3","4","5","6","7","8","9","10","11","12","13"),
      tran_cure    = c("1","2","4","5","7","8","10","12"),
      pi_init      = NULL,
      model0       = NULL,
      start        = list(fit_cox = result_elf$fit_cox,
                          fit_alpha = result_elf$fit_logit),
      epsilon_tol  = 1e-4,
      max_iter     = 200,
      t_cure       = NA
  ) 
    },
    error = function(e) {
      message("Bootstrap ", i, " failed: ", e$message)
      # Record the failure
      failed_runs <<- rbind(failed_runs, 
                            data.frame(iteration = i, error = e$message,
                                       stringsAsFactors = FALSE))
      return(NULL)  
    })
 
 if (is.null(boot_result)) {
   next   # skip the rest of the loop body
 }
 
 # -- Only for successful runs --
 results_boot[[i - N.boot1 + 1]] <- pool_est(
   alpha_est = boot_result$logit_coef,
   beta_est  = boot_result$cox_coef,
   tran_state = c("12","13","15","16","24","25","26",
                  "34","35","36","45","46","56"),
   no.iter  = boot_result$iterations
 )
 
 # Remove NULL entries (failed runs) before combining
 successful_results <- Filter(Negate(is.null), results_boot)
 bootstrap_results <- do.call(rbind, successful_results)
 
  write.csv(bootstrap_results, file = paste0("ebmt results\\Bootstrap_results.csv"), row.names = TRUE)
  if (nrow(failed_runs) > 0) {
    write.csv(failed_runs,
              file = paste0("ebmt results\\Bootstrap_failures.csv"),
              row.names = FALSE)
    cat("Number of failed runs:", nrow(failed_runs), "\n")

}

}

## calculate standard errors from bootstrap samples
se<-apply(bootstrap_results,2,sd)


