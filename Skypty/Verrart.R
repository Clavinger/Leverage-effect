#skrypt do basic stochastic volatility
library(pomp)
library(quantmod)
library(beepr)
library(doParallel)
library(doSNOW)
#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")

#ladowanie zwrotow z pliku
dane<-read.csv.zoo('C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane/wig_zwroty.csv', header = T,
                   sep=',')
dane=as.xts(dane)


############################################
#poczatek okresu badania 
data.poczatkowa='1996-01'
#koniec okresu badania
data.koncowa='2016-12'
############################################


zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep="")]



Breto.SVL_statenames <- c("H","G","Y_state")
Breto.SVL_rp_names <- c("mu_h","phi","sigma_eta","ksi","nu","theta")
#Breto.SVL_ivp_names <- c("G_0")
Breto.SVL_paramnames <- c(Breto.SVL_rp_names)
Breto.SVL_covarnames <- "covaryt"


rproc1 <- "
double beta,omega,omega2;
omega = rnorm(0,sigma_eta  * sqrt(1-G*G));
omega2 = rnorm(0, 1);
G =fmax(fmin(0.99,G),-0.99)+((2*ksi-nu)-nu*fmax(fmin(0.99,G),-0.99))+theta*sqrt((1+fmax(fmin(0.99,G),-0.99))*(1-fmax(fmin(0.99,G),-0.99)))*omega2;
beta = Y_state * sigma_eta;
H = mu_h*(1 - phi) + phi*H + beta * G  * exp(-H/2) + omega;
"
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"

rproc2.filt <- "
Y_state = covaryt;
"
Breto.SVL_rproc.sim <- paste(rproc1,rproc2.sim)
Breto.SVL_rproc.filt <- paste(rproc1,rproc2.filt)


Breto.SVL_initializer <- "
G = runif(-1,0);
H =rnorm(mu_h,sigma_eta/sqrt((1-phi*phi)));
Y_state = rnorm( 0,exp(H/2) );
"
Breto.SVL_rmeasure <- "
y=Y_state;
"

Breto.SVL_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"

Breto.SVL_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
Ttheta = log(theta);
"

Breto.SVL_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
Ttheta = exp(theta);
"

Breto.SVL.filt <- pomp(data=data.frame(y=as.vector(zwroty),
                                       time=1:length(zwroty)),
                       statenames=Breto.SVL_statenames,
                       paramnames=Breto.SVL_paramnames,
                       covarnames=Breto.SVL_covarnames,
                       times="time",
                       t0=0,
                       covar=data.frame(covaryt=c(0,as.vector(zwroty)),
                                        time=0:length(zwroty)),
                       tcovar="time",
                       rmeasure=Csnippet(Breto.SVL_rmeasure),
                       dmeasure=Csnippet(Breto.SVL_dmeasure),
                       rprocess=discrete.time.sim(step.fun=Csnippet(Breto.SVL_rproc.filt),delta.t=1),
                       initializer=Csnippet(Breto.SVL_initializer),
                       toEstimationScale=Csnippet(Breto.SVL_toEstimationScale), 
                       fromEstimationScale=Csnippet(Breto.SVL_fromEstimationScale)
)



##testowa wersja parametr?w
params_test <- c(
  mu_h = -0.21,       
  phi = .98,     
  sigma_eta = .1550,
  ksi=0.035,
  nu=0.1,
  theta=.01
)

pf1 <- pfilter(Breto.SVL.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)


par(mfrow=c(3,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
rho=as.xts(pf1@filter.traj[2,1,2:(dim(pf1@filter.traj)[3])],order.by=index(zwroty))
plot(rho, minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))


run_level <- 1 
Breto.SVL_Np <-          c(100,1e3,1e3)
Breto.SVL_Nmif <-        c(10, 100,150)
Breto.SVL_Nreps_eval <-  c(4,  10, 10)
Breto.SVL_Nreps_local <- c(10, 20, 10)
Breto.SVL_Nreps_global <-c(10, 20, 10)

Breto.SVL.list<-list(Breto.SVL_Np ,Breto.SVL_Nmif,Breto.SVL_Nreps_eval,
                     Breto.SVL_Nreps_local,Breto.SVL_Nreps_global )

Breto.SVL_rw.sd_rp <- 0.02
Breto.SVL_rw.sd_ivp <- 0.1
Breto.SVL_cooling.fraction.50 <- 0.5


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)




t.if.Breto.SVL <- system.time({
  if.Breto.SVL <- foreach(i=1:Breto.SVL.list[[5]][run_level] ,
                          .packages='pomp', .combine=c,.export = "Breto.SVL.list", 
                          .options.multicore=list(set.seed=TRUE)) %dopar% try(
                            pomp::mif2(Breto.SVL.filt,start=params_test,Np=Breto.SVL.list[[1]][run_level] , Nmif=Breto.SVL.list[[2]][run_level] ,cooling.type="geometric",
                                       cooling.fraction.50=Breto.SVL_cooling.fraction.50,
                                       transform=TRUE,
                                       rw.sd = rw.sd(
                                         mu_h      = Breto.SVL_rw.sd_rp,
                                         phi       = Breto.SVL_rw.sd_rp,
                                         sigma_eta = Breto.SVL_rw.sd_rp,
                                         ksi = 0.01*Breto.SVL_rw.sd_rp,
                                         nu = 0.01*Breto.SVL_rw.sd_rp,
                                         theta= 0.01*Breto.SVL_rw.sd_rp
                                       )
                            )
                            
                          )
  
  L.if.Breto.SVL <- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list", 
                            .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                              logmeanexp(
                                replicate(Breto.SVL.list[[3]][run_level] ,
                                          logLik(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL[[i]]),Np=Breto.SVL.list[[1]][run_level]  ))
                                ),
                                se=TRUE)
                            )
  
  H.if.Breto.SVL<- foreach(i=1:Breto.SVL.list[[5]][run_level] ,.packages='pomp', .export = "Breto.SVL.list", 
                           .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                             exp(pfilter(Breto.SVL.filt,params=coef(if.Breto.SVL[[i]]),Np=Breto.SVL.list[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                           )
})


stopCluster(cl)
beep(2)
plot(if.Breto.SVL)
