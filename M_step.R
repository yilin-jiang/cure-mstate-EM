library(data.table)

calculate_weight_mstate <- function(data_longlong, id, tran_noncure, tran_cure, tran_correction,t_cure){
  
  dt <- as.data.table(data_longlong)
  
  # likelihood products
  dt[, Prob_D_g1 := prod(lik[trans %in% tran_cure]), by = id]
  dt[, Prob_D_g0 := prod(lik[trans %in% tran_noncure]), by = id]
  
  # correction 1: if event occurred in correction transitions → cannot be cured
  dt[, has_event := any(status == 1 & trans %in% tran_correction), by = id]
  dt[has_event == TRUE, Prob_D_g1 := 0]
  
  # Correction 2: if no event in tran_correction AND time to censoring > t_cure → Prob_D_g0 = 0
  if (!is.na(t_cure)) {
    # Find the maximum time for each subject (censoring time)
    dt[, max_time := max(Tstop, na.rm = TRUE), by = id]
    
    # Check conditions for each subject
    dt[, censored_after_t_cure := max_time > t_cure]
    
    # Apply Prob_D_g0 = 0 for subjects meeting both conditions
    dt[has_event == FALSE & censored_after_t_cure == TRUE, Prob_D_g0 := 0]
   
     # Clean up temporary columns (only created inside if statement)
    dt[, `:=`(max_time = NULL, censored_after_t_cure = NULL)]
  }
  
  
  # expectation
  dt[, Expect_G_D := (Prob_D_g1 * pi) /
       (Prob_D_g1 * pi + Prob_D_g0 * (1 - pi))]
  
  # weights
  dt[trans %in% tran_cure, weight := Expect_G_D]
  dt[trans %in% tran_noncure, weight := 1 - Expect_G_D]
  
  # Clean up temporary columns
  dt[, has_event := NULL]
  
  return(as.data.frame(dt))
}