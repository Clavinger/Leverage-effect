#skrypt do basic stochastic volatility z u?cyiem pakietu RStan

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
data.poczatkowa='2001-01'
#koniec okresu badania
data.koncowa='2016-12'
############################################

#zwroty do badania
zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep="")]
zworty_dat=list(T=length(zwroty),
                y=as.vector(zwroty))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Skypty")
fit <- stan(file = 'bsv.stan', data =zworty_dat , 
            iter = 1000, warmup=100, chains = 2, cores = 2)

setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")
#save(fit, file="bsv_stan.rda")
#load(fit, file="bsv_stan.rda")

print(fit)
plot(fit)

str(extract(fit))
hist(extract(fit)$mu)
stan_plot(fit,pars=c("mu","phi","sigma"))
stan_trace(fit,pars=c("mu","phi","sigma"),nrow = 3)
stan_scat(fit,pars=c("mu","phi"))
stan_scat(fit,pars=c("mu","sigma"))
stan_scat(fit,pars=c("phi","sigma"))
stan_hist(fit,pars=c("mu","phi","sigma"),nrow = 3)


