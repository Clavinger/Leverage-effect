#skrypt do basic stochastic volatility z u¿cyiem pakietu stochvol

library(stochvol)
library(quantmod)
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

#tworzenie ³anchow MCMC
res <- svsample(as.vector(zwroty), priormu = c(0, 100), priorphi = c(20, 1.1),
                priorsigma = 0.1)

#podsumowanie
summary(res, showlatent = FALSE)

#rysunek procesu zmniennosci
volplot(res, dates = index(zwroty))
volplot(res, forecast = 100, dates = index(zwroty))

#diagnostyka
par(mfrow = c(3, 1))
paratraceplot(res)
par(mfrow = c(1, 1))
plot(res, showobs = FALSE)

#reszty
myresid <- resid(res)
plot(myresid, zwroty)
   
#rysuki diagnostyczne
plotCumu(res$para)
plotAuto(res$para)
plotDens(res$para)
plotQuant(res$para)
plotSplom(res$para)
