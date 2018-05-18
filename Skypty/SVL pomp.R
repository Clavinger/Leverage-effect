#skrypt do levarage stochastic volatility
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



####nazwy
svl_statenames <- c("H","Y_state")
svl_rp_names <- c("mu_h","phi","sigma_eta","rho")
#svl_ivp_names <- c("H_0")
#svl_paramnames <- c(svl_rp_names,svl_ivp_names)
svl_paramnames <- c(svl_rp_names)
svl_covarnames <- "covaryt"



rproc1 <- "
double beta,omega;
omega = rnorm(0,sigma_eta * sqrt(1-rho*rho));
beta = Y_state * sigma_eta ;
H = mu_h*(1 - phi) + phi*H + beta * rho * exp(-H/2) + omega;
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
svl_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
svl_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
#H=H_0
svl_initializer <- "
H = rnorm(mu_h,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) );
"
###????
svl_rmeasure <- "
y=Y_state;
"

####rozk?ad warunkowy zmiennej Y
svl_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
svl_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
Trho = 0.5*log((rho+1)/(1-rho));
"

svl_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
Trho= (exp(2*rho)-1)/(exp(2*rho)+1);
"


####wypelnianie modelu danymi
svl.filt <- pomp(data=data.frame(y=as.vector(zwroty),
                                 time=1:length(zwroty)),
                 statenames=svl_statenames,
                 paramnames=svl_paramnames,
                 covarnames=svl_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0,as.vector(zwroty)),
                                  time=0:length(zwroty)),
                 tcovar="time",
                 rmeasure=Csnippet(svl_rmeasure),
                 dmeasure=Csnippet(svl_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(svl_rproc.filt),delta.t=1),
                 initializer=Csnippet(svl_initializer),
                 toEstimationScale=Csnippet(svl_toEstimationScale), 
                 fromEstimationScale=Csnippet(svl_fromEstimationScale)
)


plot(svl.filt)




##testowa wersja parametr?w
params_test <- c(
  mu_h = -0.21,       
  phi = .98,     
  rho=-0.5,
  sigma_eta = .1550
)


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 3

#liczba czasteczek
svl_Np <-          c(1000,1e3,1e3)
svl_Nmif <-        c(10, 50,150)
svl_Nreps_eval <-  c(4,  10,  20)
svl_Nreps_local <- c(4, 10, 20)
svl_Nreps_global <-c(4, 10, 10)

svllist<-list(svl_Np ,svl_Nmif,svl_Nreps_eval,
              svl_Nreps_local,svl_Nreps_global )

#parametry do metody mif2
svl_rw.sd_rp <- 0.02
svl_rw.sd_ivp <- 0.1
svl_cooling.fraction.50 <- 0.5


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)





t.if.svl <- system.time({
  if.svl <- foreach(i=1:svllist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "svllist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(svl.filt,start=params_test,Np=svllist[[1]][run_level] , Nmif=svllist[[2]][run_level] ,cooling.type="geometric",
                                 cooling.fraction.50=svl_cooling.fraction.50,
                                 transform=TRUE,
                                 rw.sd = rw.sd(
                                   mu_h      = svl_rw.sd_rp,
                                   phi       = svl_rw.sd_rp,
                                   rho       = svl_rw.sd_rp,
                                   sigma_eta = svl_rw.sd_rp
                                 )
                      )
                      
                    )
  
  L.if.svl <- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist", 
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        logmeanexp(
                          replicate(svllist[[3]][run_level] ,
                                    logLik(pfilter(svl.filt,params=coef(if.svl[[i]]),Np=svllist[[1]][run_level]  ))
                          ),
                          se=TRUE)
                      )
  
  H.if.svl<- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist", 
                     .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                       exp(pfilter(svl.filt,params=coef(if.svl[[i]]),Np=svllist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                     )
})


stopCluster(cl)
plot(if.svl)
beep(1)
#save(if.svl,L.if.svl, H.if.svl, file="svl_if_eval.rda")

r.if.svl <- data.frame(logLik=L.if.svl[,1],logLik_se=L.if.svl[,2],t(sapply(if.svl,coef)))
summary(r.if.svl$logLik,digits=5)
r.if.svl[which.max(r.if.svl$logLik),]
pairs(~logLik+mu_h+phi+sigma_eta,data=r.if.svl)


############################################################################
############################################################################
############################################################################
#rysunki procesu zmiennosci

params_nowe<- c(
  mu_h        = as.numeric(coef(if.svl[which.max(r.if.svl$logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef(if.svl[which.max(r.if.svl$logLik)])[1,'phi']), 
  rho         = as.numeric(coef(if.svl[which.max(r.if.svl$logLik)])[1,'rho']),
  sigma_eta   = as.numeric(coef(if.svl[which.max(r.if.svl$logLik)])[1,'sigma_eta'])
)

pf1 <- pfilter(svl.filt,params=params_nowe,
               Np=svllist[[1]][run_level],filter.traj=T)

par(mfrow=c(2,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(2,1))
############################################################################
############################################################################
############################################################################



##--------Likelihood maximization using randomized starting values--------

svl_box <- rbind(
  sigma_eta=c(0.01,0.3),
  phi    = c(0.9,0.99),
  mu_h = c(-1,1),
  rho= c(-1,0)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


t.svl.box <- system.time({
  if.svl.box <- foreach(i=1:svllist[[5]][run_level],.packages='pomp', .export = "svllist",.combine=c,
                        .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(
      svl.filt, Np=svllist[[1]][run_level] , Nmif=svllist[[2]][run_level] ,
      start=apply(svl_box,1,function(x) runif(1,x[1],x[2])),cooling.type="geometric",
      cooling.fraction.50=svl_cooling.fraction.50,
      transform=TRUE,
      rw.sd = rw.sd(
        mu_h      = svl_rw.sd_rp,
        phi       = svl_rw.sd_rp,
        rho       = svl_rw.sd_rp,
        sigma_eta = svl_rw.sd_rp
      )
    )
  
  L.svl.box <- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist",.combine=rbind,
                       .options.multicore=list(set.seed=TRUE)) %dopar% {
                         set.seed(87932)
                         logmeanexp(
                           replicate(svllist[[3]][run_level] ,
                                     logLik(pfilter(svl.filt,params=coef(if.svl.box[[i]]),Np=svllist[[1]][run_level] ))
                           ), 
                           se=TRUE)
                       }
  
  H.svl.box<- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist", 
                      .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        exp(pfilter(svl.filt,params=coef(if.svl.box[[i]]),Np=svllist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                      )
})



stopCluster(cl)
beep(2)
plot(if.svl.box )


#save(if.svl.box,L.svl.box,H.svl.box, file="svl_box_eval.rda")


r.box <- data.frame(logLik=L.svl.box [,1],logLik_se=L.svl.box [,2],t(sapply(if.svl.box,coef)))
summary(r.box$logLik,digits=5)
round(r.box [which.max(r.box $logLik),],4)
apply(r.box,MARGIN = 2,mean)
############################################################################
############################################################################
############################################################################
#rysunki procesu zmiennosci

params_nowe2<- c(
  mu_h        = as.numeric(coef( if.svl.box  [which.max(r.box $logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef( if.svl.box [which.max(r.box $logLik)])[1,'phi']),  
  rho         = as.numeric(coef( if.svl.box [which.max(r.box $logLik)])[1,'rho']), 
  sigma_eta   = as.numeric(coef( if.svl.box [which.max(r.box $logLik)])[1,'sigma_eta']) 
)

pf2 <- pfilter(svl.filt,params=params_nowe2,
               Np=svllist[[1]][run_level],filter.traj=T)

par(mfrow=c(2,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf2@filter.traj[1,1,2:(dim(pf2@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))

############################################################################
############################################################################
############################################################################


############################################################################
############################################################################
############################################################################
#profile funkcji wiarygodnosci


params_nowe2

#parametr mu_h
xx1<-seq(from=params_nowe2['mu_h']-1,to=params_nowe2['mu_h']+1,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.svl.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "svllist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logmeanexp(
                        replicate(svllist[[3]][run_level] ,
                                  logLik(pfilter(svl.filt,params=c(mu_h=xx1[i], 
                                                                   phi=as.numeric(params_nowe2['phi']),
                                                                   sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                   rho=as.numeric(params_nowe2['rho'])),
                                                 Np=svllist[[1]][run_level] ))
                        ), 
                        se=FALSE)
                    }

stopCluster(cl)
beep(1)

plot(xx1, L.svl.log, type='l',xlab=expression(mu[h]),ylab="logLik")
points(xx1, L.svl.log)
points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red')


#parametr phi
xx2<-seq(from=0.95,to=0.999,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.svl.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "svllist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logmeanexp(
                         replicate(svllist[[3]][run_level] ,
                                   logLik(pfilter(svl.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                    phi=xx2[i],
                                                                    sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                    rho=as.numeric(params_nowe2['rho'])),
                                                  Np=svllist[[1]][run_level] ))
                         ), 
                         se=FALSE)
                     }

stopCluster(cl)
beep(2)

plot(xx2, L.svl.log2, type='l',xlab=expression(phi),ylab="logLik")
points(xx2, L.svl.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red')



#parametr sigma_eta
xx3<-seq(from=params_nowe2['sigma_eta']-.1,to=params_nowe2['sigma_eta']+.1,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.svl.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "svllist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logmeanexp(
                         replicate(svllist[[3]][run_level] ,
                                   logLik(pfilter(svl.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                    phi=as.numeric(params_nowe2['phi']),
                                                                    sigma_eta=xx3[i],
                                                                    rho=as.numeric(params_nowe2['rho'])),
                                                  Np=svllist[[1]][run_level] ))
                         ), 
                         se=FALSE)
                     }

stopCluster(cl)
beep(1)

plot(xx3, L.svl.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
points(xx3, L.svl.log3)
points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red')

#parametr rho
xx4<-seq(from=params_nowe2['rho']-.3,to=params_nowe2['rho']+.3,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.svl.log4<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "svllist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logmeanexp(
                         replicate(svllist[[3]][run_level] ,
                                   logLik(pfilter(svl.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), 
                                                                    phi=as.numeric(params_nowe2['phi']),
                                                                    sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                                    rho=xx4[i]),
                                                  Np=svllist[[1]][run_level] ))
                         ), 
                         se=FALSE)
                     }

stopCluster(cl)
beep(2)

plot(xx4, L.svl.log4, type='l',xlab=expression(rho),ylab="logLik")
points(xx4, L.svl.log4)
points(r.box[,'rho'], r.box[,'logLik'] ,col='red')
