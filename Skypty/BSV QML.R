require("FKF")


#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")

#parametry
mu = -10
phi = 0.98
sigma = 0.2

load(file ='dane do przykladu2')

y=log(sim1.sim@data^2)
n=length(y)

OUss <- function(mu,phi,sigma){
  Tt <- matrix(phi,ncol=1)
  Zt <- matrix(1,ncol=1)
  ct <- matrix(-1.27,ncol=1)
  dt <- matrix(mu*(1-phi), nrow = 1)
  GGt<- matrix(pi^2/2,nrow = 1,ncol=1)
  HHt<- matrix(sigma^2,nrow=1,ncol=1)
  a0 <-  as.vector(matrix(mu+1.27))
  P0 <- matrix(1,nrow=1,ncol=1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
}

KF.log <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(y, nrow=1,ncol=n))
  return(-ans$logLik)
}



KF <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(y, nrow=1,ncol=n))
  return(ans$att[1,])
}




KF.log(c(mu,phi,sigma))

plot(exp(sim1.sim@states[1,]/2),type='l')
lines(exp(KF(c(mu,phi,sigma))/2+1.27),col='red')

plot(sim1.sim@states[1,],type='l')
lines((KF(c(mu,phi,sigma))+1.27),col='red')


KF.opt<-optim(c(mu,phi,sigma),KF.log, 
               method="L-BFGS-B",hessian = T,lower=c(-Inf,-1,0),upper=c(Inf,1,Inf))

round(KF.opt$par,4)
(KF.opt$value+log(2*pi)*n/2)/3498.405


