library(stochvol)
library(plotMCMC)
library(forecast)
library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
#library(nimble)

#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")

#dlugosc szeregu czasowego
n=1000
#parametry
mu = -0.5
phi = 0.98
sigma = 0.2

## Simulate a highly persistent SV process of length 500
#sim <- svsim(n, mu = mu, phi = phi, sigma = sigma)

#save(sim,file ='dane do przykladu')
load(file ='dane do przykladu2')
plot(sim1.sim)
y=sim1.sim@data

#tworzenie ?anchow MCMC
start_time <- Sys.time()
res <- svsample(y,  draws = 10000, burnin = 1000,
                priormu = c(0, 100), priorphi = c(20, 1.1),
                priorsigma = 0.1)
end_time <- Sys.time()
end_time - start_time


#podsumowanie
summary(res, showlatent = FALSE)
round(mean(res$ para[,1] ),4)
round(mean(res$ para[,2] ),4)
round(mean(res$ para[,3] ),4)

#diagnostyka
par(mfrow = c(3, 1))
paratraceplot(res)
par(mfrow = c(1, 1))
plot(res, showobs = FALSE)

#rysuki diagnostyczne
plotCumu(res$para)
plotAuto(res$para)
plotDens(res$para)
plotQuant(res$para)
plotSplom(res$para)

################################################################################
################################################################################
################################################################################
################################################################################
#POMP
################################################################################
################################################################################



zwroty= sim1.sim@data

####nazwy
bsv_statenames <- c("H","Y_state")
bsv_rp_names <- c("mu_h","phi","sigma_eta")
#bsv_ivp_names <- c("H_0")
#bsv_paramnames <- c(bsv_rp_names,bsv_ivp_names)
bsv_paramnames <- c(bsv_rp_names)
bsv_covarnames <- "covaryt"



rproc1 <- "
double omega;
omega = rnorm(0,sigma_eta );
H = mu_h*(1 - phi) + phi*H + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"
###do wypelniania danych
rproc2.filt <- "
Y_state = covaryt;
"

###symulacja modelu SVL
bsv_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
bsv_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
#H=H_0
bsv_initializer <- "
H = rnorm(mu_h  ,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) );
"
###????
bsv_rmeasure <- "
y=Y_state;
"

####rozk?ad warunkowy zmiennej Y
bsv_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
bsv_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
"

bsv_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
"


####wypelnianie modelu danymi
bsv.filt <- pomp(data=data.frame(y=as.vector(zwroty),
                                 time=1:length(zwroty)),
                 statenames=bsv_statenames,
                 paramnames=bsv_paramnames,
                 covarnames=bsv_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0,as.vector(zwroty)),
                                  time=0:length(zwroty)),
                 tcovar="time",
                 rmeasure=Csnippet(bsv_rmeasure),
                 dmeasure=Csnippet(bsv_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(bsv_rproc.filt),delta.t=1),
                 initializer=Csnippet(bsv_initializer),
                 toEstimationScale=Csnippet(bsv_toEstimationScale), 
                 fromEstimationScale=Csnippet(bsv_fromEstimationScale)
)


plot(bsv.filt)
##testowa wersja parametr?w
params_test <- c(
  mu_h = mu,       
  phi = phi,     
  sigma_eta =sigma
)


#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(bsv.filt,params=params_test,
               Np=1000,filter.traj=T)
par(mfrow=c(1,1))
plot(exp(sim1.sim@states[1,]/2),type='l',ylab=expression(sigma),xlab="time")
lines(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l',col='red',lty=2)
legend(x='topleft',legend = c('true','filtered'),
       col=c('black','red'),lty=c(1,2))
plot(pf1)
pf1$loglik

###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 3

#liczba czasteczek
bsv_Np <-          c(1000,1e3,1e3)
bsv_Nmif <-        c(10, 50,150)
bsv_Nreps_eval <-  c(4,  10,  20)
bsv_Nreps_local <- c(4, 10, 10)
bsv_Nreps_global <-c(4, 10, 10)

bsvlist<-list(bsv_Np ,bsv_Nmif,bsv_Nreps_eval,
              bsv_Nreps_local,bsv_Nreps_global )

#parametry do metody mif2
bsv_rw.sd_rp <- 0.02
bsv_rw.sd_ivp <- 0.1
bsv_cooling.fraction.50 <- 0.5




detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)




start_time <- Sys.time()
t.if.bsv <- system.time({
  if.bsv <- foreach(i=1:bsvlist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "bsvlist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(bsv.filt,start=params_test,Np=bsvlist[[1]][run_level] , Nmif=bsvlist[[2]][run_level] ,cooling.type="geometric",
                                 cooling.fraction.50=bsv_cooling.fraction.50,
                                 transform=TRUE,
                                 rw.sd = rw.sd(
                                   mu_h      = bsv_rw.sd_rp,
                                   phi       = bsv_rw.sd_rp,
                                   sigma_eta = bsv_rw.sd_rp
                                 )
                      )
                      
                    )
  
  L.if.bsv <- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist", 
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        logmeanexp(
                          replicate(bsvlist[[3]][run_level] ,
                                    logLik(pfilter(bsv.filt,params=coef(if.bsv[[i]]),Np=bsvlist[[1]][run_level]  ))
                          ),
                          se=TRUE)
                      )
  
  H.if.bsv<- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist", 
                     .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                       exp(pfilter(bsv.filt,params=coef(if.bsv[[i]]),Np=bsvlist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                     )
})


stopCluster(cl)
end_time <- Sys.time()
end_time - start_time
beep(2)
plot(if.bsv)


r.if.bsv <- data.frame(logLik=L.if.bsv[,1],logLik_se=L.if.bsv[,2],t(sapply(if.bsv,coef)))
summary(r.if.bsv$logLik,digits=5)
r.if.bsv[which.max(r.if.bsv$logLik),]
pairs(~logLik+mu_h+phi+sigma_eta,data=r.if.bsv)



bsv_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu_h = c(-1,01)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.bsv.box <- system.time({
  if.bsv.box <- foreach(i=1:bsvlist[[5]][run_level],.packages='pomp', .export = "bsvlist",.combine=c,
                        .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(bsv.filt,Np=bsvlist[[1]][run_level] , Nmif=bsvlist[[2]][run_level] ,
      start=apply(bsv_box,1,function(x) runif(1,x[1],x[2])),cooling.type="geometric", cooling.fraction.50=bsv_cooling.fraction.50,
                  transform=TRUE,
                  rw.sd = rw.sd(
                    mu_h      = bsv_rw.sd_rp,
                    phi       = bsv_rw.sd_rp,
                    sigma_eta = bsv_rw.sd_rp  )
    )
  L.bsv.box <- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                       .options.multicore=list(set.seed=TRUE)) %dopar% {
                         set.seed(87932)
                         logmeanexp(
                           replicate(bsvlist[[3]][run_level] ,
                                     logLik(pfilter(bsv.filt,params=coef(if.bsv.box[[i]]),Np=bsvlist[[1]][run_level] ))
                           ), 
                           se=TRUE)
                       }
  
  H.bsv.box<- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist", 
                      .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        exp(pfilter(bsv.filt,params=coef(if.bsv.box[[i]]),Np=bsvlist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                      )
})



stopCluster(cl)
end_time <- Sys.time()
end_time - start_time
beep(2)
plot(if.bsv.box )


#save(if.bsv.box,L.bsv.box,H.bsv.box, file="bsv_box_eval.rda")
#load( file="bsv_box_eval.rda")

r.box <- data.frame(logLik=L.bsv.box [,1],logLik_se=L.bsv.box [,2],t(sapply(if.bsv.box,coef)))
round(apply(as.matrix(r.box),MARGIN=2, FUN=mean ),4)
round(r.box [which.max(r.box $logLik),],4)




############################################################################
############################################################################
############################################################################
#rysunki procesu zmiennosci

params_nowe2<- c(
  mu_h        = as.numeric(coef( if.bsv.box  [which.max(r.box $logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef( if.bsv.box [which.max(r.box $logLik)])[1,'phi']),    
  sigma_eta   = as.numeric(coef( if.bsv.box [which.max(r.box $logLik)])[1,'sigma_eta']) 
)

############################################################################
############################################################################
############################################################################
#wykres filtracji vs symulaja POMP
pf1 <- pfilter(bsv.filt,params=params_nowe2,
               Np=bsvlist[[1]][run_level],filter.traj=T)
pf2 <- pfilter(bsv.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l',col='red',lty=2,
     ylab=expression(sigma),xlab="time")
#lines(100*exp(pf2@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l',col='blue',lty=3)
plot(exp(sim1.sim@states[1,]/2),type='l',ylab=expression(sigma),xlab="time")
lines(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),type='l',col='red',lty=2)
legend(x='topleft',legend = c('true','filtered'),
       col=c('black','red'),lty=c(1,2))

legend(x='topright',legend = c('true','filtered'),
       col=c('black','red'),lty=c(1,2))

accuracy(100*exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),100*sim$vol)
accuracy(100*exp(pf2@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),100*sim$vol)

############################################################################
############################################################################
############################################################################
#profile funkcji wiarygodnosci


params_nowe2

#parametr mu_h
xx1<-seq(from=-1.5,to=0.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logLik(pfilter(bsv.filt,params=c(mu_h=xx1[i], phi=as.numeric(params_nowe2['phi']),sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                                 Np=bsvlist[[1]][run_level] ))
                    }

stopCluster(cl)
beep(2)

plot(xx1, L.bsv.log, type='l',xlab=expression(mu[h]),ylab="logLik")
#abline(v=params_nowe2['mu_h'],lty=2)
#points(xx1, L.bsv.log)
points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red',pch=16)
p=loess(L.bsv.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)

#parametr phi
xx2<-seq(from=0.95,to=0.99,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                                   logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), phi=xx2[i],sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                                  Np=bsvlist[[1]][run_level] ))
                     }
stopCluster(cl)
beep(2)

plot(xx2, L.bsv.log2, type='l',xlab=expression(phi),ylab="logLik")
#points(xx2, L.bsv.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red',pch=16)
p2=loess(L.bsv.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)


#parametr sigma_eta
xx3<-seq(from=params_nowe2['sigma_eta']-.1,to=params_nowe2['sigma_eta']+.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                                   logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), phi=as.numeric(params_nowe2['phi']),sigma_eta=xx3[i]),
                                                  Np=bsvlist[[1]][run_level] ))
                     }

stopCluster(cl)
beep(2)

plot(xx3, L.bsv.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
#points(xx3, L.bsv.log3)
points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red',pch=16)
p3=loess(L.bsv.log3~xx3,span=0.5)
lines(xx3,p3$fitted,col='blue',lwd=2)



par(mfrow=c(1,3))
plot(xx1, L.bsv.log, xlab=expression(mu[h]),ylab="logLik")
#abline(v=params_nowe2['mu_h'],lty=2)
#points(xx1, L.bsv.log)
points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red',pch=16)
p=loess(L.bsv.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)

plot(xx2, L.bsv.log2,  xlab=expression(phi),ylab="logLik")
#points(xx2, L.bsv.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red',pch=16)
p2=loess(L.bsv.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)

plot(xx3, L.bsv.log3,  xlab=expression(sigma[eta]),ylab="logLik")
#points(xx3, L.bsv.log3)
points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red',pch=16)
p3=loess(L.bsv.log3~xx3,span=0.5)
lines(xx3,p3$fitted,col='blue',lwd=2)

par(mfrow=c(1,1))

params_nowe3<- c(
  mu_h        = -0.4749,
  phi         =  0.9692,
  sigma_eta   =  0.2833

)

params_nowe4<- c(
  mu_h        = -0.4326,
  phi         =  0.975,
  sigma_eta   =  0.2417
)

params_nowe5<- c(
  mu_h        = -0.0927,
  phi         =  0.9599,
  sigma_eta   =  0.1249

)


set.seed(87932)
logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), phi=as.numeric(params_nowe2['phi']),sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
               Np=bsvlist[[1]][run_level] ))

set.seed(87932)
logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe3['mu_h']), phi=as.numeric(params_nowe3['phi']),sigma_eta=as.numeric(params_nowe3['sigma_eta'])),
Np=bsvlist[[1]][run_level] ))

set.seed(87932)
logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe4['mu_h']), phi=as.numeric(params_nowe4['phi']),sigma_eta=as.numeric(params_nowe4['sigma_eta'])),
               Np=bsvlist[[1]][run_level] ))
set.seed(87932)
logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe5['mu_h']), phi=as.numeric(params_nowe5['phi']),sigma_eta=as.numeric(params_nowe5['sigma_eta'])),
               Np=bsvlist[[1]][run_level] ))






#parametr mu_h
xx1<-seq(from=-1.5,to=0.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logLik(pfilter(bsv.filt,params=c(mu_h=xx1[i], phi=as.numeric(params_nowe2['phi']),sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                     Np=bsvlist[[1]][run_level] ))
                    }

stopCluster(cl)
beep(2)




############################################################################
############################################################################
############################################################################
#Liu and West 
############################################################################
############################################################################

stochVCode <- nimbleCode({
  
  
  x[1] ~ dnorm(mu*(1-phi) + phi * x0, var = 1/sigmaSquaredInv)
  y[1] ~ dnorm(0, var =  exp(x[1]))
  for(t in 2:T){
    x[t] ~ dnorm(mu*(1-phi) + phi * x[t-1],  var = 1/sigmaSquaredInv)
    y[t] ~ dnorm(0, var = exp(x[t]))
  }
  
  x0 ~ dnorm(0, var = 1/sigmaSquaredInv)
  phi <- 2 * phiStar - 1
  phiStar ~ dbeta(18, 1)
  sigmaSquaredInv ~ dgamma(5, 20)
  mu ~ dnorm(0, var = 100)      
})

stochVolModel <- nimbleModel(code = stochVCode, name ='stochVol',
                             constants = list(T = length(zwroty)),  data = list(y = as.vector(zwroty)),
                             inits = list(mu = mu,  phiStar = (phi+1)/2,
                                          sigmaSquaredInv = 1/sigma^2))


sigmaSquaredSamples <- 1/sqrt((as.matrix(CstochVolLiuWestFilter$mvEWSamples,'sigmaSquaredInv')))/10
plot(density(sigmaSquaredSamples))
muSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'mu')
plot(density(muSamples))
phiSamples <- (as.matrix(CstochVolLiuWestFilter$mvEWSamples,'phiStar')+1)/2
plot(density(phiSamples))
round(mean(sigmaSquaredSamples),4)
round(mean(muSamples),4)
round(mean(phiSamples),4)

params_nowe3<- c(
  mu_h        = mean(muSamples) ,  
  phi         = mean(phiSamples),   
  sigma_eta   = mean(sigmaSquaredSamples)
)

logLik(pfilter(bsv.filt,params=params_nowe3 ,Np= 1000))

############################################################################
############################################################################
############################################################################
#PMCM
############################################################################
############################################################################
############################################################################
hyperparams <- list(min = c(-20,0.9,0), max = c(0,1,1) )
bsv.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}


pmcmc1 <-   pmcmc(pomp(bsv.filt, dprior = bsv.dprior), start = params_nowe2,
                  Nmcmc = 2000, Np = 100, max.fail = Inf,
                  proposal = mvn.diag.rw(c(mu_h = 0.01, phi = 0.01, sigma_eta = 0.01)))

continue( pmcmc1 ,Nmcmc=5000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
plot(  pmcmc1 )

plot(  pmcmc1 )
coef(pmcmc1 )
logLik(pmcmc1 )
beep(2)


