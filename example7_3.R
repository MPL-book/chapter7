## Application to a dementia dataset

# The following code is closely based on the vigenette 
# which can be found at
# <cran.r-project.org/web/packages/survivalMPLdc/vignettes/README.html?>

library(survivalMPLdc)

# set up PRIME data set
data(PRIME)
names(PRIME)<-c("Time", "Institutionalized", "Drop Out",
                "Age", "Gender", "High Education", "Alzheimer Disease",
                "Baseline CDR", "Baseline MMSE", "Baseline SMAF",
                "Baseline ZBI", "Baseline NPI", "Benzodiazepine",
                "Antipsychotic", "Living Alone", "3-month Change MMSE",
                "3-month Change SMAF", "3-month Change NPI"
)

#time observed, event and dependent censoring indicators
surv<-as.matrix(PRIME[,1:3])

#covariates 
cova<-as.matrix(PRIME[, -c(1:3)])

# the number covariates
n<-dim(PRIME)[1]
p<-dim(PRIME)[2]-3

# MPL under independent censoring
binCount <- 30
control <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = binCount, tie = 'Yes',
                                tau = 0, copula = 'independent',
                                pent = 'penalty_mspl', smpart = 'REML',
                                penc = 'penalty_mspl', smparc = 'REML',
                                mid = 1, asy = 1, ac = 1, cv = 1,
                                cat.smpar = 'No')


coxMPLests_tau0 <- coxph_mpl_dc(surv=surv,
                                cova=cova, control=control, )

MPL_beta_tau0<-coef(object = coxMPLests_tau0,
                    parameter = "beta")
MPL_phi_tau0<-coef(object = coxMPLests_tau0,
                   parameter = "phi")
mpl_h0t_tau0 <- coxMPLests_tau0$mpl_h0t
mpl_h0Ti_tau0 <- approx( surv[,1], mpl_h0t_tau0,
                         xout = seq(0, max(surv[,1]), 0.01),
                         method="constant", rule = 2, ties = mean)$y

# MPL estimate under tau=0.2
control <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = binCount, tie = 'Yes', tau = 0.2,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 'REML', penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
coxMPLests_tau0.2 <- coxph_mpl_dc(surv=surv,
                                  cova=cova, control=control, )

MPL_beta_tau0.2<-coef(object =
                        coxMPLests_tau0.2, parameter = "beta")
MPL_phi_tau0.2<-coef(object = coxMPLests_tau0.2,
                     parameter = "phi")

mpl_h0t_tau0.2 <- coxMPLests_tau0.2$mpl_h0t
mpl_h0Ti_tau0.2 <- approx( surv[,1], mpl_h0t_tau0.2,
                           xout = seq(0, max(surv[,1]), 0.01),
                           method="constant", rule = 2, ties = mean)$y

# MPL estimate under tau=0.5
control <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = binCount, tie = 'Yes', tau = 0.5,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 'REML',penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
coxMPLests_tau0.5 <- coxph_mpl_dc(surv=surv,
                                  cova=cova, control=control, )
mpl_beta_phi_zp_tau0.5 <-
  coxMPLests_tau0.5$mpl_beta_phi_zp
MPL_beta_tau0.5<-coef(object =
                        coxMPLests_tau0.5, parameter = "beta")
MPL_phi_tau0.5<-coef(object = coxMPLests_tau0.5,
                     parameter = "phi")
mpl_h0t_tau0.5 <- coxMPLests_tau0.5$mpl_h0t
mpl_h0Ti_tau0.5 <- approx( surv[,1],
                           mpl_h0t_tau0.5, xout = seq(0, max(surv[,1]), 0.01),
                           method="constant", rule = 2,
                           ties = mean)$y

# MPL estimate under tau=0.8
control <- coxph_mpl_dc.control(ordSp = 4,
                                binCount = binCount, tie = 'Yes', tau = 0.8,
                                copula = 'frank', pent = 'penalty_mspl',
                                smpart = 'REML', penc = 'penalty_mspl',
                                smparc = 'REML', mid = 1, asy = 1, ac = 1,
                                cv = 1, cat.smpar = 'No' )
coxMPLests_tau0.8 <- coxph_mpl_dc(surv=surv,
                                  cova=cova, control=control, )
mpl_beta_phi_zp_tau0.8 <-
  coxMPLests_tau0.8$mpl_beta_phi_zp
MPL_beta_tau0.8<-coef(object =
                        coxMPLests_tau0.8, parameter = "beta")
MPL_phi_tau0.8<-coef(object = coxMPLests_tau0.8,
                     parameter = "phi")
mpl_h0t_tau0.8 <- coxMPLests_tau0.8$mpl_h0t
mpl_h0Ti_tau0.8 <- approx( surv[,1],
                           mpl_h0t_tau0.8, xout = seq(0, max(surv[,1]), 0.01),
                           method="constant", rule = 2, ties = mean)$y