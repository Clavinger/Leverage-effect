#skrypt do basic stochastic volatility z uzyciem RStan

library(rstan)
library(quantmod)
library(rstudioapi)
library(plotMCMC)
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

#zwroty do badania
zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep="")]
zworty_dat=list(T=length(zwroty),
                y=as.vector(zwroty))

rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Skypty")
fit <- stan(file = 'svl.stan', data =zworty_dat , 
            iter = 10000, warmup=1000, chains = 1)

setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")
#save(fit, file="svl_stan.rda")
#load(fit, file="svl_stan.rda")
print(fit)
show(fit)
stan_plot(fit,pars=c("mu","phi","sigma","rho"))
stan_trace(fit,pars=c("mu","phi","sigma","rho"),nrow = 4)
stan_hist(fit,pars=c("mu","phi","sigma","rho"),nrow = 4)
pairs(fit)

