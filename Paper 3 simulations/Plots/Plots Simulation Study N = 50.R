#The following code is to combine and arrange the output from the simulation study
#according to type of design (accuracy, location, no effect) for the purpose of 
#creating plots and investigating the performance.

B0I <- c(3,1,0,-1,-3)
B0II <- c(3,1,0,-1,-3)
B1I <- c(2,1,0.5,0,-0.5,-1,-2)
B1II <- c(2,1,0.5,0,-0.5,-1,-2)

#Create matrices indicating if a design is a location/any effect design
LocArray <- array(NA, dim = c(500, 10, 800))
LocIArray <- array(NA, dim = c(500, 10, 256))
AccArray <- array(NA, dim = c(500, 10, 144))
NoEffArray <- array(NA, dim = c(500, 10, 25))

#Create output matrices for performance measures for  all designs
LocationI <- matrix(NA, 256, 56)
Location <- matrix(NA, 800, 56)
Accuracy <- matrix(NA, 144, 56)
NoEff <- matrix(NA, 25, 56)

colnames(Location) <- c(rep(c("B0I", "B1I", "B0II", "B1II", "ax", "ac", "Bc", "mBc", "BcmX"), 6), "AI", "NEI")
colnames(LocationI) <- c(rep(c("B0I", "B1I", "B0II", "B1II", "ax", "ac", "Bc", "mBc", "BcmX"), 6), "AI", "NEI")
colnames(Accuracy) <- c(rep(c("B0I", "B1I", "B0II", "B1II", "ax", "ac", "Bc", "mBc", "BcmX"), 6), "AI", "NEI")
colnames(NoEff) <- c(rep(c("B0I", "B1I", "B0II", "B1II", "ax", "ac", "Bc", "mBc", "BcmX"), 6), "AI", "NEI")


AcmeasLoc <- matrix(NA, 800, 6)
AcmeasLocI <- matrix(NA, 256, 6)
AcmeasA <- matrix(NA, 144, 6)
AcmeasNoEff <- matrix(NA, 25, 6)

colnames(AcmeasLoc) <- c(rep("SSDO",6))
colnames(AcmeasLocI) <- c(rep("SSDO",6))
colnames(AcmeasA) <- c(rep("SSDO",6))
colnames(AcmeasNoEff) <- c(rep("SSDO",6))

#Create matrices for the proportions of zero effects for Bc, SAM (BcmX) and AS (mBc)
PZeroEffBcLoc <- rep(NA, 800)
PZeroEffBcLocI <- rep(NA, 256)
PZeroEffBcA <- rep(NA, 144)
PZeroEffmBcLoc <- rep(NA, 800)
PZeroEffmBcLocI <- rep(NA, 256)
PZeroEffmBcA <- rep(NA, 144)
PZeroEffBcmXLoc <- rep(NA, 800)
PZeroEffBcmXLocI <- rep(NA, 256)
PZeroEffBcmXA <- rep(NA, 144)



#Set counts
count <- 0
countAccuracy <- 0
countNoEff <- 0
countColumn <- 0
countLocation <- 0
countLocationI <- 0
IsAccuracy <- 0
IsNoEff <- 0

for(i in 1:7){
  for(j in 1:7){
    for(k in 1:5){
      for(z in 1:5){
        
        count <- count + 1
        
        Sumres <- read.table(file=paste("Paper 3 simulations/ResultsN=50/result", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".txt", sep=""), sep="\t")
        #load results for whether HPD SSDO contains 0
        load(file = paste("Paper 3 simulations/ResultsN=50/Accuracyindicator", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        #load results for whether HPDs of either B1I or B1II do not contain 0
        load(file = paste("Paper 3 simulations/ResultsN=50/NoEffindicator", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        #load results for whether HPDs of Bc, AS (mBC) and SAM (BxmX) contain 0
        load(file = paste("Paper 3 simulations/ResultsN=50/ZeroEffBc", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        load(file = paste("Paper 3 simulations/ResultsN=50/ZeroEffmBc", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        load(file = paste("Paper 3 simulations/ResultsN=50/ZeroEffBcmX", B0I[k], ",", B1I[i], ",", B0II[z], ",", B1II[j], ".Rdata", sep=""))
        
        SDisOrigin <- abs(sqrt((B0I[k]+B1I[i]*Sumres[5,1])^2 + (B0II[z]+B1II[j]*Sumres[5,1])^2))
        
        if(B1I[i] == B1II[j] & B1I[i] == 0){
          IsNoEff <- 1
        }else{
          IsNoEff <- 0
          if(B0II[z] == 0 | B0I[k] == 0 | B1II[j] == 0 | B1I[i] == 0){
            if(B0II[z] == 0 & B0I[k] == 0){
              IsAccuracy <- 1
            } else if((B0I[k] == B1I[i]) & B0I[k] == 0){
              IsAccuracy <- 1
            } else if((B0II[z] == B1II[j]) & B0II[z] ==0){
              IsAccuracy <- 1
            } else if((B1I[i]/B1II[j] - B0I[k]/B0II[z]) == 0){
              IsAccuracy <- 1
            }else{
              IsAccuracy <- 0
            }
          }else{
            if((B1I[i]/B1II[j] - B0I[k]/B0II[z]) == 0){
              IsAccuracy <- 1
            }else{
              IsAccuracy <- 0
            }
          }
        }
        
        if(IsAccuracy == 1 & IsNoEff == 0){
          
          countAccuracy <- countAccuracy + 1
          
          AccArray[,1,countAccuracy] <- Accuracyindicator
          AccArray[,2,countAccuracy] <- NoEffindicator
          AccArray[,3,countAccuracy] <- ZeroEffBc
          AccArray[,4,countAccuracy] <- ZeroEffmBc
          AccArray[,5,countAccuracy] <- ZeroEffBcmX
          
          for(ii in 1:500){
            
            if(ZeroEffBc[ii] == ZeroEffmBc[ii]){
              AccArray[ii,6,countAccuracy] = 1
            }else{
              AccArray[ii,6,countAccuracy] = 0
            }
            
            if(ZeroEffBcmX[ii] == ZeroEffmBc[ii]){
              AccArray[ii,7,countAccuracy] = 1
            }else{
              AccArray[ii,7,countAccuracy] = 0
            }
            
            if(ZeroEffBc[ii] == ZeroEffBcmX[ii]){
              AccArray[ii,8,countAccuracy] = 1
            }else{
              AccArray[ii,8,countAccuracy] = 0
            }
            
            if(NoEffindicator[ii] == ZeroEffBc[ii]){
              AccArray[ii,9,countAccuracy] = 1
            }else{
              AccArray[ii,9,countAccuracy] = 0
            }
            
            if(Accuracyindicator[ii] == ZeroEffBc[ii]){
              AccArray[ii,10,countAccuracy] = 0
            }else{
              AccArray[ii,10,countAccuracy] = 1
            }
          }
          
          for(s in 1:6){
            for(p in 1:9){
              
              countColumn <- countColumn + 1
              Accuracy[countAccuracy, countColumn] <- Sumres[p,s]
              
            }
            AcmeasA[countAccuracy, s] <- Sumres[10,s]
          }
          
          Accuracy[countAccuracy, 55] <- sum(Accuracyindicator)/500
          Accuracy[countAccuracy, 56] <- sum(NoEffindicator)/500
          PZeroEffBcA[countAccuracy] <- sum(ZeroEffBc)/500
          PZeroEffmBcA[countAccuracy] <- sum(ZeroEffmBc)/500
          PZeroEffBcmXA[countAccuracy] <- sum(ZeroEffBcmX)/500
          countColumn <- 0
          
        }
        
        
        if(IsAccuracy == 0 & IsNoEff == 0){
          
          if(SDisOrigin > 1){
            
            countLocation <- countLocation + 1
            
            LocArray[,1, countLocation] <- Accuracyindicator
            LocArray[,2, countLocation] <- NoEffindicator
            LocArray[,3, countLocation] <- ZeroEffBc
            LocArray[,4, countLocation] <- ZeroEffmBc
            LocArray[,5, countLocation] <- ZeroEffBcmX
            
            for(ii in 1:500){
              if(ZeroEffBc[ii] == ZeroEffmBc[ii]){
                LocArray[ii,6,countLocation] = 1
              }else{
                LocArray[ii,6,countLocation] = 0
              }
              
              if(ZeroEffBcmX[ii] == ZeroEffmBc[ii]){
                LocArray[ii,7,countLocation] = 1
              }else{
                LocArray[ii,7,countLocation] = 0
              }
              
              if(ZeroEffBc[ii] == ZeroEffBcmX[ii]){
                LocArray[ii,8,countLocation] = 1
              }else{
                LocArray[ii,8,countLocation] = 0
              }
              
              if(NoEffindicator[ii] == ZeroEffBc[ii]){
                LocArray[ii,9,countLocation] = 1
              }else{
                LocArray[ii,9,countLocation] = 0
              }
              
              if(Accuracyindicator[ii] == ZeroEffBc[ii]){
                LocArray[ii,10,countLocation] = 0
              }else{
                LocArray[ii,10,countLocation] = 1
              }
            }
            
            for(s in 1:6){
              for(p in 1:9){
                
                countColumn <- countColumn + 1
                Location[countLocation, countColumn] <- Sumres[p,s]
                
              }
              AcmeasLoc[countLocation, s] <- Sumres[10,s]
            }
            
            Location[countLocation, 55] <- sum(Accuracyindicator)/500
            Location[countLocation, 56] <- sum(NoEffindicator)/500
            PZeroEffBcLoc[countLocation] <- sum(ZeroEffBc)/500
            PZeroEffmBcLoc[countLocation] <- sum(ZeroEffmBc)/500
            PZeroEffBcmXLoc[countLocation] <- sum(ZeroEffBcmX)/500
            countColumn <- 0
            
          }else if(SDisOrigin <= 1){
            
            countLocationI <- countLocationI + 1
            
            LocIArray[,1, countLocationI] <- Accuracyindicator
            LocIArray[,2, countLocationI] <- NoEffindicator
            LocIArray[,3, countLocationI] <- ZeroEffBc
            LocIArray[,4, countLocationI] <- ZeroEffmBc
            LocIArray[,5, countLocationI] <- ZeroEffBcmX
            
            for(ii in 1:500){
              if(ZeroEffBc[ii] == ZeroEffmBc[ii]){
                LocIArray[ii,6,countLocationI] = 1
              }else{
                LocIArray[ii,6,countLocationI] = 0
              }
              
              if(ZeroEffBcmX[ii] == ZeroEffmBc[ii]){
                LocIArray[ii,7,countLocationI] = 1
              }else{
                LocIArray[ii,7,countLocationI] = 0
              }
              
              if(ZeroEffBc[ii] == ZeroEffBcmX[ii]){
                LocIArray[ii,8,countLocationI] = 1
              }else{
                LocIArray[ii,8,countLocationI] = 0
              }
              
              if(NoEffindicator[ii] == ZeroEffBc[ii]){
                LocIArray[ii,9,countLocationI] = 1
              }else{
                LocIArray[ii,9,countLocationI] = 0
              }
              
              if(Accuracyindicator[ii] == ZeroEffBc[ii]){
                LocIArray[ii,10,countLocationI] = 0
              }else{
                LocIArray[ii,10,countLocationI] = 1
              }
            }
            
            for(s in 1:6){
              for(p in 1:9){
                
                countColumn <- countColumn + 1
                LocationI[countLocationI, countColumn] <- Sumres[p,s]
                
              }
              AcmeasLocI[countLocationI, s] <- Sumres[10,s]
            }
            
            LocationI[countLocationI, 55] <- sum(Accuracyindicator)/500
            LocationI[countLocationI, 56] <- sum(NoEffindicator)/500
            PZeroEffBcLocI[countLocationI] <- sum(ZeroEffBc)/500
            PZeroEffmBcLocI[countLocationI] <- sum(ZeroEffmBc)/500
            PZeroEffBcmXLocI[countLocationI] <- sum(ZeroEffBcmX)/500
            countColumn <- 0
          }
        }
        
        
        if(IsNoEff == 1){
          
          countNoEff <- countNoEff + 1
          
          NoEffArray[,1,countNoEff] <- Accuracyindicator
          NoEffArray[,2,countNoEff] <- NoEffindicator
          NoEffArray[,3,countNoEff] <- ZeroEffBc
          NoEffArray[,4,countNoEff] <- ZeroEffmBc
          NoEffArray[,5,countNoEff] <- ZeroEffBcmX
          
          for(ii in 1:500){
            if(ZeroEffBc[ii] == ZeroEffmBc[ii]){
              NoEffArray[ii,6,countNoEff] = 1
            }else{
              NoEffArray[ii,6,countNoEff] = 0
            }
            
            if(ZeroEffBcmX[ii] == ZeroEffmBc[ii]){
              NoEffArray[ii,7,countNoEff] = 1
            }else{
              NoEffArray[ii,7,countNoEff] = 0
            }
            
            if(ZeroEffBc[ii] == ZeroEffBcmX[ii]){
              NoEffArray[ii,8,countNoEff] = 1
            }else{
              NoEffArray[ii,8,countNoEff] = 0
            }
            
            if(NoEffindicator[ii] == ZeroEffBc[ii]){
              NoEffArray[ii,9,countNoEff] = 1
            }else{
              NoEffArray[ii,9,countNoEff] = 0
            }
            
            if(Accuracyindicator[ii] == ZeroEffBc[ii]){
              NoEffArray[ii,10,countNoEff] = 0
            }else{
              NoEffArray[ii,10,countNoEff] = 1
            }
          }
          
          for(s in 1:6){
            for(p in 1:9){
              
              countColumn <- countColumn + 1
              NoEff[countNoEff, countColumn] <- Sumres[p,s]
              
            }
            AcmeasNoEff[countNoEff, s] <- Sumres[10,s]
          }
          
          NoEff[countNoEff, 55] <- sum(Accuracyindicator)/500
          NoEff[countNoEff, 56] <- sum(NoEffindicator)/500
          countColumn <- 0
          
        }
        
        
        
      }
    }
  }
}


library(ggplot2)
library(RColorBrewer)
library(extrafont)
font_install("fontcm")
loadfonts()
loadfonts(device = "win")
par(family = "LM Roman 10")

#Data for crosstabulations 5-7
Tabulate <- function(x, y){
  X <- as.factor(x)
  levels(X) <- c(0,1)
  Y <- as.factor(y)
  levels(Y) <- c(0,1)
  return(table(X,Y))
}

TabAcc <- (sapply((1:144), function(w){Tabulate(AccArray[,3,w], AccArray[,1,w])})/500)*100
TabLoc <- (sapply((1:800), function(w){Tabulate(LocArray[,3,w], LocArray[,1,w])})/500)*100
TabLocI <- (sapply((1:256), function(w){Tabulate(LocIArray[,3,w], LocIArray[,1,w])})/500)*100

TabBAcc <- (sapply((1:144), function(w){c(sum(AccArray[,2,w]), 500-sum(AccArray[,2,w]))})/500)*100
TabBLoc <- (sapply((1:800), function(w){c(sum(LocArray[,2,w]), 500-sum(LocArray[,2,w]))})/500)*100
TabBLocI <- (sapply((1:256), function(w){c(sum(LocIArray[,2,w]), 500-sum(LocIArray[,2,w]))})/500)*100

rowMeans(TabAcc)
rowMeans(TabLoc)
rowMeans(TabLocI)

rowMeans(TabBAcc)
rowMeans(TabBLoc)
rowMeans(TabBLocI)

library(matrixStats)
rowSds(TabAcc)
rowSds(TabLoc)
rowSds(TabLocI)
rowSds(TabBAcc)
rowSds(TabBLoc)
rowSds(TabBLocI)




################Plots for Relative Bias, Coverage and AIW#######################
pdf("Paper 3 simulations/Plots/FigureBiasAcmeasBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,25]) / abs(Location[,7]),
                   abs(LocationI[,25]) / abs(LocationI[,7]),
                   abs(Accuracy[,7] - Accuracy[,25]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
       scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
     theme(legend.position = "none",
           text = element_text(family="CM Roman", size = 24), 
           axis.title.x  = element_text(size = 24),
           axis.title.y  = element_text(size = 24),
           plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("Relative Bias\n")  
dev.off()

png("Paper 3 simulations/Plots/FigureBiasAcmeasBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,25]) / abs(Location[,7]),
                   abs(LocationI[,25]) / abs(LocationI[,7]),
                   abs(Accuracy[,7] - Accuracy[,25]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
  theme(legend.position = "none",
        text = element_text(family="LM Roman 10", size = 16), 
        axis.title.x  = element_text(size = 16),
        axis.title.y  = element_text(size = 16),
        plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Relative Bias\n")  
dev.off()

pdf("Paper 3 simulations/Plots/FigureBiasAcmeasmBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,26]) / abs(Location[,8]),
                   abs(LocationI[,26]) / abs(LocationI[,8]),
                   abs(Accuracy[,8] - Accuracy[,26]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
       scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
     theme(legend.position = "none",
           text = element_text(family="CM Roman", size = 24), 
           axis.title.x  = element_text(size = 24),
           axis.ticks.y  = element_blank(),
           axis.text.y = element_blank(),
           axis.title.y  = element_blank(),
           plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("Relative Bias\n")
dev.off()

png("Paper 3 simulations/Plots/FigureBiasAcmeasmBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,26]) / abs(Location[,8]),
                   abs(LocationI[,26]) / abs(LocationI[,8]),
                   abs(Accuracy[,8] - Accuracy[,26]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
  theme(legend.position = "none",
        text = element_text(family="LM Roman 10", size = 16), 
        axis.title.x  = element_text(size = 16),
        axis.ticks.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
        plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Relative Bias\n")
dev.off()

pdf("Paper 3 simulations/Plots/FigureBiasAcmeasBcmX.pdf", family = "CM Roman", pointsize = 24)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,27]) / abs(Location[,9]),
                   abs(LocationI[,27]) / abs(LocationI[,9]),
                   abs(Accuracy[,9] - Accuracy[,27]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
       scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
     theme(legend.position = "none",
           text = element_text(family="CM Roman", size = 24), 
           axis.title.x  = element_text(size = 24),
           axis.ticks.y  = element_blank(),
           axis.text.y = element_blank(),
           axis.title.y  = element_blank(),
           plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Relative Bias\n") 
dev.off()

png("Paper 3 simulations/Plots/FigureBiasAcmeasBcmX.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD  <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Bias          <- c(abs(Location[,27]) / abs(Location[,9]),
                   abs(LocationI[,27]) / abs(LocationI[,9]),
                   abs(Accuracy[,9] - Accuracy[,27]))
Type          <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d   <- data.frame(Bias, PopulationSD, Type)
dg  <- qplot(PopulationSD, Bias, colour = Type, data = d) + ylim(0, 0.8) +
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))+
  theme(legend.position = "none",
        text = element_text(family="LM Roman 10", size = 16), 
        axis.title.x  = element_text(size = 16),
        axis.ticks.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y  = element_blank(),
        plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Relative Bias\n") 
dev.off()

pdf("Paper 3 simulations/Plots/FigureCoverageAcmeasBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,52], LocationI[,52], Accuracy[,52])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.title.y  = element_text(size = 24),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureCoverageAcmeasBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,52], LocationI[,52], Accuracy[,52])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.title.y  = element_text(size = 16),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

pdf("Paper 3 simulations/Plots/FigureCoverageAcmeasmBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,53], LocationI[,53], Accuracy[,53])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureCoverageAcmeasmBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,53], LocationI[,53], Accuracy[,53])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()


pdf("Paper 3 simulations/Plots/FigureCoverageAcmeasBcmX.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,54], LocationI[,54], Accuracy[,54])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureCoverageAcmeasBcmX.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
Cov          <- c(Location[,54], LocationI[,54], Accuracy[,54])
Type         <- as.factor(c(rep("Location",1056), rep("Accuracy",144)))

d <- data.frame(Cov, PopulationSD, Type) 
dg <- qplot(PopulationSD,Cov, colour = Type, data=d) + ylim(0.85,1) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("Coverage\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

pdf("Paper 3 simulations/Plots/FigureAIWAcmeasBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,43] - Location[,34]) / abs(Location[,16]),
                      abs(LocationI[,43] - LocationI[,34]) / abs(LocationI[,16]),
                      abs(Accuracy[,43] - Accuracy[,34]) / abs(Accuracy[,16])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.title.y  = element_text(size = 24),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) + 
     xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureAIWAcmeasBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,43] - Location[,34]) / abs(Location[,16]),
                      abs(LocationI[,43] - LocationI[,34]) / abs(LocationI[,16]),
                      abs(Accuracy[,43] - Accuracy[,34]) / abs(Accuracy[,16])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.title.y  = element_text(size = 16),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) + 
  xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

pdf("Paper 3 simulations/Plots/FigureAIWAcmeasmBc.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,44] - Location[,35]) / abs(Location[,17]),
                      abs(LocationI[,44] - LocationI[,35]) / abs(LocationI[,17]),
                      abs(Accuracy[,44] - Accuracy[,35]) / abs(Accuracy[,17])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureAIWAcmeasmBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,44] - Location[,35]) / abs(Location[,17]),
                      abs(LocationI[,44] - LocationI[,35]) / abs(LocationI[,17]),
                      abs(Accuracy[,44] - Accuracy[,35]) / abs(Accuracy[,17])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

pdf("Paper 3 simulations/Plots/FigureAIWAcmeasBcmX.pdf", family = "CM Roman", pointsize = 24)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,45] - Location[,36]) / abs(Location[,18]),
                      abs(LocationI[,45] - LocationI[,36]) / abs(LocationI[,18]),
                      abs(Accuracy[,45] - Accuracy[,36]) / abs(Accuracy[,18])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
      scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="CM Roman", size = 24), 
                             axis.title.x  = element_text(size = 24),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
     xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

png("Paper 3 simulations/Plots/FigureAIWAcmeasBcmX.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
PopulationSD <- c(AcmeasLoc[,1], AcmeasLocI[,1], AcmeasA[,1])
AIW          <- log(c(abs(Location[,45] - Location[,36]) / abs(Location[,18]),
                      abs(LocationI[,45] - LocationI[,36]) / abs(LocationI[,18]),
                      abs(Accuracy[,45] - Accuracy[,36]) / abs(Accuracy[,18])))
Type         <- as.factor(c(rep("Location", 1056), rep("Accuracy", 144)))

d <- data.frame(AIW, PopulationSD, Type)
dg <- qplot(PopulationSD, AIW, colour = Type, data = d) + ylim(-3,12) + 
  scale_color_manual(breaks = c("Location", "Accuracy"), values=c("black", "grey"))

dg + theme_classic() + theme(legend.position = "none",
                             text=element_text(family="LM Roman 10", size = 16), 
                             axis.title.x  = element_text(size = 16),
                             axis.ticks.y  = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y  = element_blank(),
                             plot.margin = unit(c(0,0,0.5,0.5), "cm")) +
  xlab("\nTrue SSDO") + ylab("log(AIW)\n") +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        panel.grid.major.y=element_line(color="grey", size = 0.5))
dev.off()

###############Histograms of posterior modes for 3 datasets#####################
require(Rcpp)
sourceCpp("Paper 3 simulations/SimulationCode.cpp")

set.seed(101)

B0I <- 3
B0II <- 1
B1I <- 2
B1II <- 0.5
N <- 50
p <- 2
mu <- 0
sd <- 1
tm <- 5000
t.lag <- 1
burn <- 0
nsim <- 500
EXBurn <- 1000

Res <- SimRegRCPP1P(N, p, c(B0I, B1I), c(B0II, B1II), mu, sd, tm, t.lag, burn, nsim, EXBurn)

Max <- c()
Mac <- c()
MBc <- c()
MmBc <- c()
MBcmX <- c()

for(i in 1:nsim){
  Max[i] <- hmode(Res$ax[1001:5000,,i], 0.1)
  Mac[i] <- hmode(Res$ac[1001:5000,,i], 0.1)
  MBc[i] <- hmode(Res$Bc[1001:5000,,i], 0.1)
  MmBc[i] <- hmode(Res$mBc[1001:5000,,i], 0.1)
  MBcmX[i] <- hmode(Res$BcmX[1001:5000,,i], 0.1)
}

hist(Max)
hist(Mac)
hist(MBc)
hist(MmBc)
hist(MBcmX)


pdf("Paper 3 simulations/Plots/FigureHist3210,5ac.pdf", family = "CM Roman", pointsize = 26)
par(mar = c(4,4,0.1,0) + 0.1)
hist(Mac%%(2*pi), main= "", xlab = expression(a[c]), ylim = c(0,150), breaks = 60)
dev.off()

png("Paper 3 simulations/Plots/FigureHist3210,5ac.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
par(mar = c(4,4,0.1,0) + 0.1)
hist(Mac%%(2*pi), main= "", xlab = expression(a[c]), ylim = c(0,150), breaks = 60)
dev.off()

pdf("Paper 3 simulations/Plots/FigureHist3210,5Bc.pdf", family = "CM Roman", pointsize = 26)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MBc, main= "", xlab = expression(b[c]), ylab = "", ylim = c(0,150), yaxt='n', breaks = 60)
dev.off()

png("Paper 3 simulations/Plots/FigureHist3210,5Bc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MBc, main= "", xlab = expression(b[c]), ylab = "", ylim = c(0,150), yaxt='n', breaks = 60)
dev.off()

pdf("Paper 3 simulations/Plots/FigureHist3210,5mBc.pdf", family = "CM Roman", pointsize = 26)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MmBc, main= "", xlab = "AS", ylim = c(0,150), breaks = 60)
dev.off()

png("Paper 3 simulations/Plots/FigureHist3210,5mBc.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MmBc, main= "", xlab = "AS", ylim = c(0,150), breaks = 60)
dev.off()

pdf("Paper 3 simulations/Plots/FigureHist3210,5BcmX.pdf", family = "CM Roman", pointsize = 26)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MBcmX, main= "", xlab = "SAM", ylab = "", ylim = c(0,150), yaxt='n', breaks = 60)
dev.off()

png("Paper 3 simulations/Plots/FigureHist3210,5BcmX.png", family = "LM Roman 10", width = 5, height = 5, units = "in", pointsize = 16, res = 1200)
par(mar = c(4,4,0.1,0) + 0.1)
hist(MBcmX, main= "", xlab = "SAM", ylab = "", ylim = c(0,150), yaxt='n', breaks = 60)
dev.off()

library(zoo)

###################Convergence Plot of Exemplary Design#########################

pdf("Paper 3 simulations/Plots/FigureConvergence3210,5.pdf", 
    height = 4, width = 8, family = "CM Roman", pointsize = 12)
data <- as.matrix(cbind(Res$ac[4500:5000,,300]%%(2*pi),Res$Bc[4500:5000,,300]))
plot(as.zoo(data), ylab = expression(a[c], b[c]), main = "", xlab = "Iteration")
dev.off()

png("Paper 3 simulations/Plots/FigureConvergence3210,5.png", 
    height = 4, width = 8, family = "LM Roman 10", units = "in", pointsize = 12, res = 1200)
data <- as.matrix(cbind(Res$ac[4500:5000,,300]%%(2*pi),Res$Bc[4500:5000,,300]))
plot(as.zoo(data), ylab = expression(a[c], b[c]), main = "", xlab = "Iteration")
dev.off()

