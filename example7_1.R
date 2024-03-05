## Generate random data with dependent censoring
library(copula)

n=2000
#define the copula
tau_ec <- 0.7
alpha<-iTau(claytonCopula(100), tau = tau_ec)
cop <- claytonCopula(alpha, dim = 2)

x<-rep(0, n) #observed time
del<-rep(0, n) #event status
eta<-rep(0, n) #(dependent) censoring status
y <- NULL
c <- NULL

#generate covariates
Z <- cbind(rbinom(n,1,0.5), runif(n,-1,1)) 

#set true parameter values
betas <- c(1,-0.5)
phis <- c(0.2, -1)

#set weibull distribution parameters
lambda <- 2
a <- 14/3
#lamda and a can be selected differently by the reader

#exp(xTb) term for event time
expbetT <- exp(Z%*%betas)
#scale of weibull baseline hazard for event
lamT<-lambda*(exp(-Z%*%betas))^(1/a) 
#scale of weibull baseilne hazard for censoring
lamCw<-lambda*(exp(-Z%*%phis))^(1/a) 

for(i in 1:n){
  S_tc <- rCopula(1, cop) #S_event and S_cens probabilities
  
  #draw random event time from weibull dist
  ti <- qweibull((1-S_tc[1,1]), a, lamT[i]) 
  #draw random censoring time from weibull dist
  ci<-qweibull((1-S_tc[1,2]), a, lamCw[i]) 
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
head(dat)

#plot
plot(y ~ c, xlab = "Censoring time", ylab = "Event time")
#fit normal Cox model (ignoring dependent censoring)
library(survival)
surv.obj <- Surv(time = x, event = del)
coxph(surv.obj ~ Z[,1] + Z[,2])
