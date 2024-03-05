## A simulation study with dependent censoring

# function for simulating data

gen_depcen <- function(n, tau){
  
  #set up copula
  alpha<-iTau(frankCopula(100), tau = tau)
  cop <- frankCopula(alpha, dim = 2)
  
  x<-rep(0, n) #observed time
  del<-rep(0, n) #event status
  eta<-rep(0, n) #(dependent) censoring status
  y <- NULL
  c <- NULL
  
  #generate covariates
  Z <- cbind(rbinom(n,1,0.5), runif(n,-10,10)) 
  
  #set true parameter values
  betas <- c(-0.5,0.1)
  phis <- c(0.3, 0.2)
  
  #set weibull distribution parameters
  lambda <- 2
  a <- 2
  
  for(i in 1:n){
    S_tc <- rCopula(1, cop) 
    #S_event and S_cens probabilities
    
    #draw random event time
    ti <- (-log(S_tc[1])/((lambda^(-a))*exp(Z[i,]%*%betas)))^(1/a)
    #draw random censoring time
    ci <- (-log(S_tc[2])/((lambda^(-a))*exp(Z[i,]%*%phis)))^(1/a)
    #save event and censoring times
    y <- c(y,ti)
    c <- c(c,ci)
    
    if(ti<ci){
      #event observed
      del[i]<-1
      x[i]<-ti
    }else if (ci<ti){
      #dependent censoring time observed
      eta[i]<-1
      x[i]<-ci
    }
  }
  dat<-cbind(x, del, eta)
  
  out <- list(dat = dat, Z = Z)
  return(out)
  
}

# function for running simulation study with n=500
ch7_save <- matrix(0, nrow = 300, ncol = 48)
library(survivalMPLdc)

for(s in 1:300){
  
  dat <- gen_depcen(n=500, tau = 0.8)
  surv.dat <- dat$dat
  Z.dat <- dat$Z
  
  #fit independent censoring
  ctrl1 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 30, tau = 0,
                                copula = 'independent', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
  coxMPLests_tau0 <- coxph_mpl_dc(surv=surv.dat,
                                  cova=Z.dat, control=ctrl1)
  
  ch7_save[s,1] <- coxMPLests_tau0$mpl_beta[1]
  ch7_save[s,2] <- coxMPLests_tau0$mpl_beta[2]
  ch7_save[s,3] <- coxMPLests_tau0$mpl_beta_sd[1]
  ch7_save[s,4] <- coxMPLests_tau0$mpl_beta_sd[2]
  
  ch7_save[s,5] <- coxMPLests_tau0$mpl_phi[1]
  ch7_save[s,6] <- coxMPLests_tau0$mpl_phi[2]
  ch7_save[s,7] <- coxMPLests_tau0$mpl_phi_sd[1]
  ch7_save[s,8] <- coxMPLests_tau0$mpl_phi_sd[2]
  
  #fit frank cop and tau = 0.8
  ctrl2 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 30, tau = 0.8,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
  coxMPLests_tau0.8_frank <- coxph_mpl_dc(surv=surv.dat,
                                          cova=Z.dat, control=ctrl2)
  
  ch7_save[s,9] <- coxMPLests_tau0.8_frank$mpl_beta[1]
  ch7_save[s,10] <- coxMPLests_tau0.8_frank$mpl_beta[2]
  ch7_save[s,11] <- coxMPLests_tau0.8_frank$mpl_beta_sd[1]
  ch7_save[s,12] <- coxMPLests_tau0.8_frank$mpl_beta_sd[2]
  
  ch7_save[s,13] <- coxMPLests_tau0.8_frank$mpl_phi[1]
  ch7_save[s,14] <- coxMPLests_tau0.8_frank$mpl_phi[2]
  ch7_save[s,15] <- coxMPLests_tau0.8_frank$mpl_phi_sd[1]
  ch7_save[s,16] <- coxMPLests_tau0.8_frank$mpl_phi_sd[2]
  
  #fit frank cop and tau = 0.4
  ctrl3 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 30, tau = 0.4,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
  coxMPLests_tau0.4_frank <- coxph_mpl_dc(surv=surv.dat,
                                          cova=Z.dat, control=ctrl3)
  
  ch7_save[s,17] <- coxMPLests_tau0.4_frank$mpl_beta[1]
  ch7_save[s,18] <- coxMPLests_tau0.4_frank$mpl_beta[2]
  ch7_save[s,19] <- coxMPLests_tau0.4_frank$mpl_beta_sd[1]
  ch7_save[s,20] <- coxMPLests_tau0.4_frank$mpl_beta_sd[2]
  
  ch7_save[s,21] <- coxMPLests_tau0.4_frank$mpl_phi[1]
  ch7_save[s,22] <- coxMPLests_tau0.4_frank$mpl_phi[2]
  ch7_save[s,23] <- coxMPLests_tau0.4_frank$mpl_phi_sd[1]
  ch7_save[s,24] <- coxMPLests_tau0.4_frank$mpl_phi_sd[2]
  
  #fit clayton cop and tau = 0.8
  ctrl4 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 30, tau = 0.8,
                                copula = 'clayton', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
  coxMPLests_tau0.8_clay <- coxph_mpl_dc(surv=surv.dat,
                                         cova=Z.dat, control=ctrl4)
  
  ch7_save[s,25] <- coxMPLests_tau0.8_clay$mpl_beta[1]
  ch7_save[s,26] <- coxMPLests_tau0.8_clay$mpl_beta[2]
  ch7_save[s,27] <- coxMPLests_tau0.8_clay$mpl_beta_sd[1]
  ch7_save[s,28] <- coxMPLests_tau0.8_clay$mpl_beta_sd[2]
  
  ch7_save[s,29] <- coxMPLests_tau0.8_clay$mpl_phi[1]
  ch7_save[s,30] <- coxMPLests_tau0.8_clay$mpl_phi[2]
  ch7_save[s,31] <- coxMPLests_tau0.8_clay$mpl_phi_sd[1]
  ch7_save[s,32] <- coxMPLests_tau0.8_clay$mpl_phi_sd[2]
  
  #fit clayton cop and tau = 0.4
  ctrl5 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 30, tau = 0.4,
                                copula = 'clayton', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
  coxMPLests_tau0.4_clay <- coxph_mpl_dc(surv=surv.dat,
                                         cova=Z.dat, control=ctrl5)
  
  ch7_save[s,33] <- coxMPLests_tau0.4_clay$mpl_beta[1]
  ch7_save[s,34] <- coxMPLests_tau0.4_clay$mpl_beta[2]
  ch7_save[s,35] <- coxMPLests_tau0.4_clay$mpl_beta_sd[1]
  ch7_save[s,36] <- coxMPLests_tau0.4_clay$mpl_beta_sd[2]
  
  ch7_save[s,37] <- coxMPLests_tau0.4_clay$mpl_phi[1]
  ch7_save[s,38] <- coxMPLests_tau0.4_clay$mpl_phi[2]
  ch7_save[s,39] <- coxMPLests_tau0.4_clay$mpl_phi_sd[1]
  ch7_save[s,40] <- coxMPLests_tau0.4_clay$mpl_phi_sd[2]
  
  #fit frank cop and tau = 0.8, maximum likelihood (no smoothing)
  ctrl6 <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = 10, tau = 0.8,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 0, penc = 'penalty_mspl',
                                smparc = 0, mid = 1, asy = 1, ac = 0,
                                cv = 0, cat.smpar = 'No' , maxit2 = 1)
  coxests_tau0.8_frank_ML <- coxph_mpl_dc(surv=surv.dat,
                                          cova=Z.dat, control=ctrl6)
  
  ch7_save[s,41] <- coxests_tau0.8_frank_ML$mpl_beta[1]
  ch7_save[s,42] <- coxests_tau0.8_frank_ML$mpl_beta[2]
  ch7_save[s,43] <- coxests_tau0.8_frank_ML$mpl_beta_sd[1]
  ch7_save[s,44] <- coxests_tau0.8_frank_ML$mpl_beta_sd[2]
  
  ch7_save[s,45] <- coxests_tau0.8_frank_ML$mpl_phi[1]
  ch7_save[s,46] <- coxests_tau0.8_frank_ML$mpl_phi[2]
  ch7_save[s,47] <- coxests_tau0.8_frank_ML$mpl_phi_sd[1]
  ch7_save[s,48] <- coxests_tau0.8_frank_ML$mpl_phi_sd[2]
  
  
}