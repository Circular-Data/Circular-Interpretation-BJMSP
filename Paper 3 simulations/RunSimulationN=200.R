#Source Simulation code and required packages for analysis
require(Rcpp)
sourceCpp("Paper 3 simulations/SimulationCode.cpp")

#Set seed
set.seed(101)

#Create vectors with population values
B0I <- c(3,1,0,-1,-3)
B0II <- c(3,1,0,-1,-3)
B1I <- c(2,1,0.5,0,-0.5,-1,-2)
B1II <- c(2,1,0.5,0,-0.5,-1,-2)
N <- 200
p <- 2
mu <- 0
sd <- 1
tm <- 5000
t.lag <- 1
burn <- 0
nsim <- 500
EXBurn <- 1000

#Start simulation
for(i in 1:7){
  for(j in 1:7){
    for(k in 1:5){
      for(z in 1:5){
        
        #Compute results for all datasets of 1 design
        Res <- SimRegRCPP1P(N, p, c(B0I[k], B1I[i]), c(B0II[z], B1II[j]), mu, sd, tm, t.lag, burn, nsim, EXBurn)
        Sumres <- Res$SumRes
        Accuracyindicator <- Res$Accuracyindicator
        NoEffindicator <- Res$NoEffindicator
        ZeroEffBc <- Res$ZeroEffBc
        ZeroEffmBc <- Res$ZeroEffmBc
        ZeroEffBcmX <- Res$ZeroEffBcmX
        #Save results and indicators
        write.table(Sumres, file=paste("Paper 3 simulations/ResultsN=200/result", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".txt", sep=""), sep="\t")
        save(Accuracyindicator, file=paste("Paper 3 simulations/ResultsN=200/Accuracyindicator", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        save(NoEffindicator, file=paste("Paper 3 simulations/ResultsN=200/NoEffindicator", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        save(ZeroEffBc, file=paste("Paper 3 simulations/ResultsN=200/ZeroEffBc", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        save(ZeroEffmBc, file=paste("Paper 3 simulations/ResultsN=200/ZeroEffmBc", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        save(ZeroEffBcmX, file=paste("Paper 3 simulations/ResultsN=200/ZeroEffBcmX", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
      }
    }
  }
}

