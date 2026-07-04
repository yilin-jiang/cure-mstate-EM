# =========================== Making scenarios for simulation ==============================#

# Find alpha_0 and censoring scales to make scenarios

# X ~ Bernoulli(0.5), V~N(0,1)##
# alpha1=0.5, alpha2=-0.2

source("generate_data.R")
# find alpha0 so that cure proportions are at 20%, 50%, 80%
set.seed(123)
X <- rbinom(1000000, 1, 0.5)
V <- rnorm(1000000, 0, 1)

cure_rate_diff <- function(alpha0, target) {
  p <- plogis(alpha0 + 0.5 * X - 0.2 * V)
  mean(p) - target
}

solve_alpha0 <- function(target) {
  uniroot(
    cure_rate_diff,
    interval = c(-5, 5),
    target = target
  )$root
}

alpha0_20 <- solve_alpha0(0.20)
alpha0_50 <- solve_alpha0(0.50)
alpha0_80 <- solve_alpha0(0.80)

# alpha0_20 <- -1.66664
# alpha0_50 <--0.2499823
# alpha0_80 <- 1.166669

# find scale_cens so that censoring rates are at 15%, 35%, 55%

find_scale_cens<-function(alpha0,
                          censor_rate){
  n_sample=10000
  covs_df<-generate_cure(n=10000,alpha0=alpha0)

  ## non-cure subjects ##
  noncure<-which(covs_df$G==0)

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

  lower=1
  upper=150
  
  iter <- 0
  cen_count <- 0
  while (abs(cen_count/ n_sample - censor_rate) > 0.01) {
    cen_count <- 0
    
    iter<-iter+1
    scale_cens<- (lower + upper) / 2
    
    for(i in noncure){
      generated_data <- icmstate::sim_weibmsm(tmat = mstate::trans.illdeath(), 
                                              shape = c(1,1,1), 
                                              scale = c(scale12[i], scale13[i], scale23[i]),
                                              n_subj =1, obs_pars = c(2, 0.5, 80),
                                              censshape = c(1,1,1),
                                              censscale=c(scale_cens,scale_cens,scale_cens),
                                              true_trajec = TRUE)
      outcome_data<-generated_data$true
      if (nrow(outcome_data)==2 & outcome_data$state[2]==1){cen_count<-cen_count+1}
      if (nrow(outcome_data)==3 & outcome_data$state[3]==2){cen_count<-cen_count+1}
    }
    
    ## cured subjects ##
    cure<-which(covs_df$G==1)
    
    lp_23_c<-as.vector(as.matrix(covs_df[,c("X","V")]) %*% c(0.4,0.1))-0.4
    lambda13_cure <- 0.05 * exp(lp_23_c)
    
    for (j in cure){
      T2 <- rexp(1,
                             rate = lambda13_cure[j])
      T_C<-rexp(1,rate=1/scale_cens)
      if (T2>T_C){cen_count<-cen_count+1}
    }
      
    if (cen_count/ n_sample < censor_rate) {
      upper <- scale_cens
    } else {
      lower <- scale_cens
    }
    print(paste0("iter ",iter))
    print(paste0("lower=",lower))
    print(paste0("upper=",upper))
    print(paste0("cens_scale=",scale_cens))
    print(paste0("censoring proportion ",cen_count/ n_sample))
  }
  return(scale_cens=scale_cens)
}

#### Make scenarios #####
cure_prop<- c("small" = 0.2, "medium" = 0.5, "high" = 0.8)
censor_rate<-c("small" = 0.15, "medium" = 0.35, "high" = 0.55)
full_factorial <- data.table::data.table(
  expand.grid(
    "cure_prop" = cure_prop,
    "censor_rate" = censor_rate
  )
) 

full_factorial <- full_factorial[, alpha0 := ifelse(cure_prop== 0.2, alpha0_20,ifelse(cure_prop==0.5,alpha0_50,alpha0_80))]
full_factorial <- full_factorial[, cen_scale:=mapply(find_scale_cens,alpha0=alpha0,censor_rate=censor_rate)]
full_factorial[, scen_num := 1:.N]
save(full_factorial, file = "scenarios.RData")
