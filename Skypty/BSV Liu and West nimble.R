library(nimble)


#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")
load(file ='dane do przykladu')

#parametry
mu = -10
phi = 0.98
sigma = 0.2

zwroty= sim$y

stochVCode <- nimbleCode({

      
      x[1] ~ dnorm(mu*(1-phi) + phi * x0, var = sigmaSquaredInv)
      y[1] ~ dnorm(0, var =  exp(x[1]))
      for(t in 2:T){
          x[t] ~ dnorm(mu*(1-phi) + phi * x[t-1],  var = sigmaSquaredInv)
          y[t] ~ dnorm(0, var = exp(x[t]))
        }
  
      x0 ~ dnorm(mu, var = sigmaSquaredInv)
      phi <- 2 * phiStar - 1
      phiStar ~ dbeta(20, 1.5)
      sigmaSquaredInv ~ dinvgamma(2.5, 0.025)
      mu ~ dnorm(0, 1/100)      
    })

stochVolModel <- nimbleModel(code = stochVCode, name ='stochVol',
                             constants = list(T = length(zwroty)),  data = list(y = zwroty),
                            inits = list(mu = mu,  phiStar = (phi+1)/2,
                                                  sigmaSquaredInv = sigma^2))
CstochVolModel <- compileNimble(stochVolModel)
stochVolLiuWestFilter <- buildLiuWestFilter(model = stochVolModel, nodes ='x', 
                                               params = c("mu","phiStar","sigmaSquaredInv"))

CstochVolLiuWestFilter <- compileNimble(stochVolLiuWestFilter,
                                           project = stochVolModel)
CstochVolLiuWestFilter$run(1000)
sigmaSquaredSamples <- sqrt((as.matrix(CstochVolLiuWestFilter$mvEWSamples,'sigmaSquaredInv')))
plot(density(sigmaSquaredSamples))
muSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'mu')
plot(density(muSamples))
phiSamples <- (as.matrix(CstochVolLiuWestFilter$mvEWSamples,'phiStar')+1)/2
plot(density(phiSamples))
round(mean(sigmaSquaredSamples),4)
round(mean(muSamples),4)
round(mean(phiSamples),4)


Rmcmc <- buildMCMC(stochVolModel)
Cmodel <- compileNimble(stochVolModel )
Cmcmc <- compileNimble(Rmcmc, project = stochVolModel )
res<-runMCMC(Cmcmc, niter = 11000, nburnin = 1000, progressBar = TRUE,samplesAsCodaMCMC = TRUE,
             nchains = 1, inits = list(mu = mu,  phiStar = (phi+1)/2,
            sigmaSquaredInv = 1/sigma^2))
str(res)
round(mean(res[,1]),4)
round(mean((2*res[,2]-1)),4)
round(mean(sqrt(res[,3])),4)
plot(density((2*res[,2]-1)))


#rysuki diagnostyczne
plotCumu(res$para)
plotAuto(res$para)
plotDens(res$para)
plotQuant(res$para)
plotSplom(res$para)
