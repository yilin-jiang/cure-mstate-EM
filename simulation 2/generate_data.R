#####Function for data generation########

library(dplyr)
library(mstate)
# function to transform data from long data format to wide
long_to_wide<-function(long_data_full,id_var="id",covars){
  
  long_data <- long_data_full %>%
    mutate(trans_label = paste0("T", from,to))  
  
  wide_data <- long_data%>%
    tidyr::pivot_wider(
      id_cols = c(!!sym(id_var), all_of(covars)),
      names_from = trans_label,
      values_from = c(Tstop, status),
      names_sep = "_"
    )
  wide_data$T2<-wide_data$Tstop_T12
  wide_data$status2<-wide_data$status_T12
  wide_data <- wide_data %>% 
    mutate(T3=pmax(Tstop_T13,Tstop_T23,na.rm=TRUE)) %>% 
    mutate(status3=pmax(status_T13,status_T23,na.rm=TRUE))%>% 
    dplyr::select(all_of(c(id_var,"T2","status2","T3","status3",covars)))
  wide_data=as.data.frame(wide_data)
  return(wide_data)
}

# function to generate covariates and cure status
generate_cure <- function(n, alpha0) {
  X <- rbinom(n, 1, 0.5)
  V <- rnorm(n,0,1)
  
  lp_cure <- alpha0 + 0.5 * X - 0.2 * V
  p   <- plogis(lp_cure)
  
  G <- rbinom(n, 1, p)
  
  data.frame(G = G, X = X, V = V)
}

# function to generate data
#Argment required: 
#n: sample size;
#alpha0:intercept value for cure regression model
#scale_cens: scale parameter for censoring distribution
generate_data <- function(n,alpha0,scale_cens){
#generate covariates and cure status
covs_df<-generate_cure(n=n,alpha0=alpha0)
covs_df$id<-1:nrow(covs_df)

## non-cure subjects ##
noncure<-covs_df$id[which(covs_df$G==0)]

# Linear predictors
lp_12 <- as.vector(as.matrix(covs_df[,c("X","V")]) %*% c(0.5,-0.3))
lp_13 <- as.vector(as.matrix(covs_df[,c("X","V")]) %*% c(0.3,0.2))
lp_23 <- as.vector(as.matrix(covs_df[,c("X","V")]) %*% c(0.4,0.1))

# Convert to scale parameters for dweibull, sim_weibmsm
scale12_0 <-1/0.1
scale13_0<-1/0.05
scale23_0<-1/0.2
scale12 <- scale12_0 * exp(-lp_12)
scale13 <- scale13_0 * exp(-lp_13)
scale23 <- scale23_0 * exp(-lp_23)

long_data_full<-data.frame("id"=integer(),"from"=integer(),"to"=integer(),
                           "trans"=integer(),"Tstart"= numeric(),"Tstop"= numeric(),
                           "time"= numeric(),"status"=integer(),
                           "X"= numeric(),"Z"= numeric())
  
  for(i in noncure){
    generated_data <- icmstate::sim_weibmsm(tmat = mstate::trans.illdeath(), 
                                            shape = c(1,1,1), 
                                            scale = c(scale12[i], scale13[i], scale23[i]),
                                            n_subj =1, obs_pars = c(2, 0.5, 80),
                                            censshape = c(1,1,1),
                                            censscale=c(scale_cens,scale_cens,scale_cens),
                                            true_trajec = TRUE)
    outcome_data<-generated_data$true
    outcome_data[,1]<-outcome_data[,1]-outcome_data[1,1] #reset entry time to state 1 as zero
    outcome_data$id<-i
    
    long_data <- data.frame(
      id     = outcome_data$id,
      from   = NA,
      to     = NA,
      trans  = NA,
      Tstart = NA,
      Tstop  = NA,
      time   = NA,
      status = NA
    )

    if (nrow(outcome_data)==2){
      long_data$from<-c(1,1)
      long_data$to<-c(2,3)
      long_data$trans<-c(1,2)
      long_data$Tstart<-outcome_data$time[1]
      long_data$Tstop<-outcome_data$time[2]
      long_data$time<-long_data$Tstop-long_data$Tstart
      if (outcome_data$state[2]==1){
        long_data$status<-c(0,0)
      }else {
        long_data$status<-c(0,1)
      }
    }else if(nrow(outcome_data)==3){
      long_data$from<-c(1,1,2)
      long_data$to<-c(2,3,3)
      long_data$trans<-c(1,2,3)
      long_data$Tstart<-c(outcome_data$time[1],outcome_data$time[1],outcome_data$time[2])
      long_data$Tstop<-c(outcome_data$time[2],outcome_data$time[2],outcome_data$time[3])
      long_data$time<-long_data$Tstop-long_data$Tstart
      if (outcome_data$state[3]==2){
        long_data$status<-c(1,0,0)
      }else {
        long_data$status<-c(1,0,1)
      }
    }
    long_data<-merge(long_data,covs_df[i,],by="id")
    long_data_full<-rbind(long_data_full,long_data)
  }

## cured subjects ##
cure<-covs_df$id[which(covs_df$G==1)]
lp_23_c<-as.vector(as.matrix(covs_df[,c("X","V")]) %*% c(0.4,0.1))-0.4

lambda13_cure <- 0.05 * exp(lp_23_c) 

for (j in cure){
  long_data <- data.frame(
    id     = rep(j,2),
    from   = rep(1,2),
    to     = c(2,3),
    trans  = c(1,2),
    Tstart = c(0,0),
    Tstop  = NA,
    time   = NA,
    status = c(0,NA)
  )
  T2 <- rexp(1, rate = lambda13_cure[j])
  T_C<-rexp(1,rate=1/scale_cens)
  long_data$Tstop<-ifelse(T2>T_C,T_C,T2)
  long_data$time<-long_data$Tstop
  long_data$status[2]<-ifelse(T2>T_C,0,1)
  
  long_data<-merge(long_data,covs_df[j,],by="id")
  long_data_full<-rbind(long_data_full,long_data)
}
long_data_full<-long_data_full[order(long_data_full$id,long_data_full$trans),]
long_data_obs<-long_data_full%>% select(-"G")
wide_data_obs<-long_to_wide(long_data_obs,covars=c("X","V"))

#convert long format data to msdata format
class(long_data_obs) <- c("msdata", "data.frame")
attr(long_data_obs, "trans") <- trans.illdeath()
return(list(long_data_full=long_data_full,long_data_obs=long_data_obs,wide_data_obs=wide_data_obs))
}

