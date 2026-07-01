library(data.table)
library(survival)
#model: a list of coxph models
#tran: transition vector that includes all possible from-to pairs, in the same order as in the model list
#data_long: data in long or extended long format with column id, from, to, time, status and covariates
calculate_hazard_dt <- function(model, tran, data_long){
  
  dt <- as.data.table(data_long)
  
  dt[, `:=`(cumHaz = NA_real_, hazard = NA_real_)]
  dt[, fromto := paste0(from, to)]
  
  for(k in seq_along(model)){
    
    mod <- model[[k]]
    
    idx <- dt[fromto == tran[k], which = TRUE]
    if(length(idx) == 0) next
    
    bh <- suppressWarnings(basehaz(mod, centered = FALSE))
    
    # Construct cumulative hazard including time 0
    bh_time <- c(0, bh$time)
    bh_cumhaz <- c(0, bh$hazard)
    
    # hazard increments (jump sizes)
    bh_jump <- diff(bh_cumhaz)
    
    # linear predictors
    lp <- predict(mod, newdata = dt[idx], type = "lp")
    risk <- exp(lp)
    
    # interval lookup
    tstop_idx  <- findInterval(dt$Tstop[idx], bh_time)
    tstart_idx <- findInterval(dt$Tstart[idx], bh_time)
    
    # cumulative hazard over interval
    cumhaz <- (bh_cumhaz[tstop_idx] - bh_cumhaz[tstart_idx]) * risk
    
    # hazard contribution at event time
    haz <- rep(1, length(idx))
    event_rows <- dt$status[idx] == 1
    
    haz[event_rows] <- bh_jump[pmax(tstop_idx[event_rows]-1,1)] * risk[event_rows]
    
    dt[idx, `:=`(
      cumHaz = cumhaz,
      hazard = haz
    )]
  }
  
  return(as.data.frame(dt))
}
