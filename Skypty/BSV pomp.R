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
data.poczatkowa='2001-01'
#koniec okresu badania
data.koncowa='2016-12'
############################################


zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep="")]


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
H = rnorm(mu_h,sigma_eta/sqrt((1-phi*phi))) ;
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
  mu_h = -0.21,       
  phi = .98,     
  sigma_eta = .1550
)


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 1

#liczba czasteczek
bsv_Np <-          c(1000,1e3,2e3)
bsv_Nmif <-        c(10, 50,200)
bsv_Nreps_eval <-  c(4,  10,  20)
bsv_Nreps_local <- c(4, 10, 20)
bsv_Nreps_global <-c(4, 10, 100)

bsvlist<-list(bsv_Np ,bsv_Nmif,bsv_Nreps_eval,
              bsv_Nreps_local,bsv_Nreps_global )

#parametry do metody mif2
bsv_rw.sd_rp <- 0.02
bsv_rw.sd_ivp <- 0.1
bsv_cooling.fraction.50 <- 0.5


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)





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
plot(if.bsv)
beep(1)
#save(if.bsv,L.if.bsv, H.if.bsv, file="bsv_if_eval.rda")

r.if.bsv <- data.frame(logLik=L.if.bsv[,1],logLik_se=L.if.bsv[,2],t(sapply(if.bsv,coef)))
summary(r.if.bsv$logLik,digits=5)
r.if.bsv[which.max(r.if.bsv$logLik),]
pairs(~logLik+mu_h+phi+sigma_eta,data=r.if.bsv)

############################################################################
############################################################################
############################################################################
#rysunki procesu zmiennosci

params_nowe<- c(
  mu_h        = as.numeric(coef(if.bsv[which.max(r.if.bsv$logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef(if.bsv[which.max(r.if.bsv$logLik)])[1,'phi']),    
  sigma_eta   = as.numeric(coef(if.bsv[which.max(r.if.bsv$logLik)])[1,'sigma_eta'])
)

pf1 <- pfilter(bsv.filt,params=params_nowe,
               Np=bsvlist[[1]][run_level],filter.traj=T)

par(mfrow=c(2,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf1@filter.traj[1,1,2:(dim(pf1@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(2,1))
############################################################################
############################################################################
############################################################################


##--------Likelihood maximization using randomized starting values--------

bsv_box <- rbind(
  sigma_eta=c(0.001,0.3),
  phi    = c(0.9,1),
  mu_h = c(-1,1)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


t.bsv.box <- system.time({
  if.bsv.box <- foreach(i=1:bsvlist[[5]][run_level],.packages='pomp', .export = "bsvlist",.combine=c,
                        .options.multicore=list(set.seed=TRUE)) %dopar%  
    pomp::mif2(
      if.bsv[[1]],
      start=apply(bsv_box,1,function(x) runif(1,x[1],x[2]))
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
plot(if.bsv.box )
beep(1)

#save(if.bsv.box,L.bsv.box,H.bsv.box, file="bsv_box_eval.rda")


r.box <- data.frame(logLik=L.bsv.box [,1],logLik_se=L.bsv.box [,2],t(sapply(if.bsv.box,coef)))
summary(r.box$logLik,digits=5)
r.box [which.max(r.box $logLik),]




############################################################################
############################################################################
############################################################################
#rysunki procesu zmiennosci

params_nowe2<- c(
  mu_h        = as.numeric(coef( if.bsv.box  [which.max(r.box $logLik)])[1,'mu_h']),    
  phi         = as.numeric(coef( if.bsv.box [which.max(r.box $logLik)])[1,'phi']),    
  sigma_eta   = as.numeric(coef( if.bsv.box [which.max(r.box $logLik)])[1,'sigma_eta']) 
)

pf2 <- pfilter(bsv.filt,params=params_nowe2,
               Np=bsvlist[[1]][run_level],filter.traj=T)

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
L.bsv.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logmeanexp(
                         replicate(bsvlist[[3]][run_level] ,
                                   logLik(pfilter(bsv.filt,params=c(mu_h=xx1[i], phi=as.numeric(params_nowe2['phi']),sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                                  Np=bsvlist[[1]][run_level] ))
                         ), 
                         se=FALSE)
                     }

stopCluster(cl)
beep(1)

plot(xx1, L.bsv.log, type='l',xlab=expression(mu[h]),ylab="logLik")
points(xx1, L.bsv.log)
points(r.box[,'mu_h'], r.box[,'logLik'] ,col='red')


#parametr phi
xx2<-seq(from=0.95,to=0.999,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logmeanexp(
                        replicate(bsvlist[[3]][run_level] ,
                                  logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), phi=xx2[i],sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                                 Np=bsvlist[[1]][run_level] ))
                        ), 
                        se=FALSE)
                    }

stopCluster(cl)
beep(1)

plot(xx2, L.bsv.log2, type='l',xlab=expression(phi),ylab="logLik")
points(xx2, L.bsv.log2)
points(r.box[,'phi'], r.box[,'logLik'] ,col='red')



#parametr sigma_eta
xx3<-seq(from=params_nowe2['sigma_eta']-.1,to=params_nowe2['sigma_eta']+.1,length.out = 10)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logmeanexp(
                         replicate(bsvlist[[3]][run_level] ,
                                   logLik(pfilter(bsv.filt,params=c(mu_h=as.numeric(params_nowe2['mu_h']), phi=as.numeric(params_nowe2['phi']),sigma_eta=xx3[i]),
                                                  Np=bsvlist[[1]][run_level] ))
                         ), 
                         se=FALSE)
                     }

stopCluster(cl)
beep(1)

plot(xx3, L.bsv.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
points(xx3, L.bsv.log3)
points(r.box[,'sigma_eta'], r.box[,'logLik'] ,col='red')



############################################################################
############################################################################
############################################################################
#PMCM
############################################################################
############################################################################
############################################################################
hyperparams <- list(min = c(-2,0,0,-1), max = c(2,1,1,1) )
bsv.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}


pmcmc1 <-   pmcmc(pomp(bsv.filt, dprior = bsv.dprior), start = params_est,
                  Nmcmc = 2000, Np = 100, max.fail = Inf,
                  proposal = mvn.diag.rw(c(mu_h = 0.01, phi = 0.01, sigma_eta = 0.01,H_0=0.01)))

continue( pmcmc1 ,Nmcmc=5000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
plot(  pmcmc1 )

plot(  pmcmc1 )
coef(pmcmc1 )
logLik(pmcmc1 )
beep(2)


