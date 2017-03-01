#Empty global environment
rm(list=ls())

#Load/source required packages/files
library(circular)
library(Rcpp)
sourceCpp('Mode Finding Algorithms/VenterMode.cpp')
sourceCpp('Mode Finding Algorithms/VenterModeCircular.cpp')
sourceCpp('Regression Sampler/RegressionRCPPextraLANmeasuresIIe.cpp')

#Read original data with absolute pointing errors
DataS1 <- read.csv("Pointing North Data/S1_Dataset.csv", sep = ";")

#####################Don't run this code twice (lengthy)!#######################
#
#Simulate data with non-absolute pointing errors based on the circular mean and 
#resultant length.
#
#repeatedly sample a vector with -1's and 1's of the same length as the sample
#size. Multiply this vector with the absolute pointing errors. Compute the 
#resultant length for this simulated non-absolute pointing error. If the 
#resultant length equals 0.43 (reported resultant length of the non-absolute 
#pointing error of the original data) stop the loop.

#set.seed(101)
#repeat{
#  toss_seq  <- sample(c(-1,1), length(DataS1$Estimate), replace = TRUE)
#  CircDat   <- as.circular(toss_seq*DataS1$Estimate, units = c("degrees"))
#  summary   <- summary(CircDat)
#  if(round(summary[8], 2) == 0.43){break}
#}

#Convert compute the simulated pointing errors and convert to degrees. Set the 
#mean to 19.57, the mean of the original data. Compute summary statistics to 
#check if procedure was done correctly (check mean and resultant length).

#newEst3   <- DataS1$Estimate*toss_seq
#newEst3_c <- as.circular(newEst3, units = c("degrees"))
#newEst3a  <- (newEst3-mean(newEst3_c)) + 19.57
#summary(as.circular(newEst3a, units = c("degrees")))

#Save the simulted data
#save(newEst3, newEst3a, toss_seq,
#     file = "Pointing North Data/toss_seq superlengthy.RData")

#Load the file with simulated data:
load(file = "Pointing North Data/toss_seq superlengthy.RData")

#Compute summary statistics of the pointing error for males and females separately
#and combined
males     <- which(DataS1$Sex..1m.2f. == 1)
females   <- which(DataS1$Sex..1m.2f. == 2)
males_c   <- as.circular(newEst3a[males], units = c("degrees"))
females_c <- as.circular(newEst3a[females], units = c("degrees"))

summary(as.circular(newEst3a, units = c("degrees")))
summary(males_c)
summary(females_c)

#Create a dataframe for analysis
North <- as.data.frame(cbind(newEst3a*(pi/180), DataS1[,1:5]))
colnames(North) <- c("theta", "Participant", "Age", "Sex", "Experience", "SBSOD")

#Obtain summary statistics for the continuous predictors
summary(North)

#Center continuous variables and recode sex to 0 = males, 1 = females.
North$Sex[North$Sex==1] <- 0
North$Sex[North$Sex==2] <- 1
North$Experience <- North$Experience - mean(North$Experience)
North$Age        <- North$Age        - mean(North$Age)
North$SBSOD      <- North$SBSOD      - mean(North$SBSOD)


#Create model matrices for analyses
X1 <- as.matrix(cbind(rep(1,length(North$theta)), North[,3:6]))
X2 <- X1

txtdata <- cbind(North$theta, X1[,-1])
colnames(txtdata) <- c("theta", "age", "sex", "experience", "sbsod")
write.table(txtdata, "Pointing North Data/DatatoSamplerPointingNorth.txt", sep = "\t", row.names = F)

#Fit a regression model and save the results
#set.seed(101)
#Results <- Regression(North$theta, X1, X2, 3000, 1, 0)
#save(Results, file = "Pointing North Data/toss_seq superlengthy Results.RData")
load(file = "Pointing North Data/toss_seq superlengthy Results.RData")

B  <- Results$B
I1 <- B[1001:3000, 1]
I2 <- B[1001:3000, 6]

#Convergence plots
plot.ts(B[,1:5])
plot.ts(B[,6:10])
plot.ts(Results$Bc)
plot.ts(Results$mBc)
plot.ts(Results$BcmX)
plot.ts(Results$SD)

#Modes and HPD for location/accuracy check
hmode(Results$SD[1001:3000, 1], 0.1)
hmode(Results$SD[1001:3000, 3], 0.1)
hmode(Results$SD[1001:3000, 4], 0.1)
hmodeci(Results$SD[1001:3000, 1], 0.95)
hmodeci(Results$SD[1001:3000, 3], 0.95)
hmodeci(Results$SD[1001:3000, 4], 0.95)

#Modes and HPD for circular effects
#Compute circular modes and HPD
modescirc <- matrix(NA, 4, 6)
HPDcirc   <- array(NA, dim = c(4, 2, 6))

colnames(modescirc) <-c("a_x", "a_c", "b_c", "b_cnew", "mbc", "bcmx")

b_cnew    <-  matrix(NA, 2000, 4)

for(i in 2:5){
  
  for(j in 1001:3000){
    
    b_cnew[j-1000, i-1]  <- atan2(B[j, 6] + B[j, i+5], B[j, 1] + B[j,i]) - atan2(B[j, 6], B[j, 1]) 
    
  }
  
  modescirc[i-1,1] <- hmode(Results$ax[1001:3000, i-1], 0.1)
  HPDcirc[i-1,,1]  <- hmodeci(Results$ax[1001:3000, i-1], 0.95)
  modescirc[i-1,2] <- hmodeC(Results$ac[1001:3000, i-1], 0.1)
  HPDcirc[i-1,,2]  <- hmodeciC(Results$ac[1001:3000, i-1], 0.95)
  
  modescirc[i-1,3] <- hmode(Results$Bc[1001:3000, i-1], 0.1)
  HPDcirc[i-1,,3]  <- hmodeci(Results$Bc[1001:3000, i-1], 0.95)
  modescirc[i-1,4] <- hmodeC(b_cnew[, i-1], 0.1)
  HPDcirc[i-1,,4]  <- hmodeciC(b_cnew[, i-1], 0.95)
  modescirc[i-1,5] <- hmode(Results$mBc[1001:3000, i-1], 0.1)
  HPDcirc[i-1,,5]  <- hmodeci(Results$mBc[1001:3000, i-1], 0.95)
  modescirc[i-1,6] <- hmode(Results$BcmX[1001:3000, i-1], 0.1)
  HPDcirc[i-1,,6]  <- hmodeci(Results$BcmX[1001:3000, i-1], 0.95)
  hist(Results$ax[1001:3000, i-1])
  plot.ts(Results$ax[1001:3000, i-1])
  hist(Results$ac[1001:3000, i-1])
  plot.ts(Results$ac[1001:3000, i-1])
  hist(Results$Bc[1001:3000, i-1])
  plot.ts(Results$Bc[1001:3000, i-1])
  hist(Results$mBc[1001:3000, i-1])
  plot.ts(Results$mBc[1001:3000, i-1])
  hist(Results$BcmX[1001:3000, i-1])
  plot.ts(Results$BcmX[1001:3000, i-1])
}


#Modes and 95% HPD for bivariate components

modeslin <- rep(NA, 10)
modesBc <- rep(NA, 4)
modesmBc <- rep(NA, 4)
modesBcmX <- rep(NA, 4)
HPDlin   <- matrix(NA, 10, 2)
HPDBc   <- matrix(NA, 4, 2)
HPDmBc   <- matrix(NA, 4, 2)
HPDBcmX   <- matrix(NA, 4, 2)

for(i in 1:5){
  
  modeslin[i]   <- hmode(B[1001:3000, i], 0.1)
  HPDlin[i,]    <- hmodeci(B[1001:3000, i], 0.95)
  modeslin[i+5] <- hmode(B[1001:3000, i+5], 0.1)
  HPDlin[i+5,]  <- hmodeci(B[1001:3000, i+5], 0.95)
}
for(i in 1:4){
  modesBc[i]    <- hmode(Results$Bc[1001:3000,i], 0.1)
  modesmBc[i]   <- hmode(Results$mBc[1001:3000,i], 0.1)
  modesBcmX[i]  <- hmode(Results$BcmX[1001:3000,i], 0.1)
  HPDBc[i,]     <- hmodeci(Results$Bc[1001:3000,i], 0.95)
  HPDmBc[i,]    <- hmodeci(Results$mBc[1001:3000,i], 0.95)
  HPDBcmX[i,]   <- hmodeci(Results$BcmX[1001:3000,i], 0.95)
}


###############################Plots############################################

library(extrafont)
font_install("fontcm")
loadfonts()

######################Posterior Histograms######################################

pdf("Pointing North Data/PostHistBrunye.pdf", height = 9, width = 6,
    family = "CM Roman", pointsize = 12)
par(mfrow = c(5,2))

hist(B[,1], breaks = 50, xlab = "", main = "Intercept")
hist(B[,6], breaks = 50, xlab = "", ylab = "", main = "")
hist(B[,2], breaks = 50, xlab = "", main = "Age")
hist(B[,7], breaks = 50, xlab = "", ylab = "", main = "")
hist(B[,3], breaks = 50, xlab = "", main = "Gender")
hist(B[,8], breaks = 50, xlab = "", ylab = "", main = "")
hist(B[,4], breaks = 50, xlab = "", main = "Experience")
hist(B[,9], breaks = 50, xlab = "", ylab = "", main = "")
hist(B[,5], breaks = 50, xlab = "", main = "SBSOD")
hist(B[,10], breaks = 50, xlab = "", ylab = "", main = "")

dev.off()

pdf("Pointing North Data/PostHistBrunyecirc.pdf", height = 9, width = 6,
    family = "CM Roman", pointsize = 12)
par(mfrow = c(4,2))

hist(Results$Bc[1001:3000,1], breaks = 50, xlab = expression(b[c]), main = "Age")
hist(Results$Bc[1001:3000,3], breaks = 50, xlab = expression(b[c]),ylab = "", main = "Experience")
hist(Results$mBc[1001:3000,1], breaks = 50, xlab = "AS",  main = "")
hist(Results$mBc[1001:3000,3], breaks = 50, xlab = "AS", ylab = "", main = "")
hist(Results$BcmX[1001:3000,1], breaks = 50, xlab = "SAM",  main = "")
hist(Results$BcmX[1001:3000,3], breaks = 50, xlab = "SAM", ylab = "", main = "")
hist(Results$SD[1001:3000,1], breaks = 50, xlab = "SSDO",  main = "")
hist(Results$SD[1001:3000,3], breaks = 50, xlab = "SSDO",  main = "")

dev.off()



###########################Convergence plots####################################
B  <- Results$B
colnames(B) <- c("Intercept", "Age", "Sex", "Experience", "SBSOD",
                 "Intercept", "Age", "Sex", "Experience", "SBSOD")

pdf("Pointing North Data/ConvergenceBrunyeB1.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)

Y <- B[,1:5]

axis2 <- function(Y, col = col, bg = bg, pch = pch, type = type, ...) {
  lines(Y, col=col)
  axis(side = 2, las = 1, at = c(round(mean(Y), 2), round(mean(Y)+sd(Y),2),round(mean(Y)-sd(Y),2)), cex.axis = 1)}

plot.ts(B[,1:5], main="", xlab = c("Iteration"), axes = F, panel = axis2)

dev.off()

pdf("Pointing North Data/ConvergenceBrunyeB2.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)

Y <- B[,6:10]

axis2 <- function(Y, col = col, bg = bg, pch = pch, type = type, ...) {
  lines(Y, col=col)
  axis(side = 2, las = 1, at = c(round(mean(Y), 2), round(mean(Y)+sd(Y),2),round(mean(Y)-sd(Y),2)), cex.axis = 1)}

plot.ts(B[,6:10], main="", xlab = c("Iteration"), axes = F, panel = axis2)

dev.off()

###############Predicted Regression Line and ConcentrationPlots#################

#SBSOD

X <- seq(-5, 5, 1)
Xreal <- seq(-2.7, 2.7, 0.1)
CompI  <- modeslin[1] + modeslin[5]*X
CompII <- modeslin[6] + modeslin[10]*X
CompIreal  <- modeslin[1] + modeslin[5]*Xreal
CompIIreal <- modeslin[6] + modeslin[10]*Xreal
axval <- hmode(Results$ax[1001:3000,4], 0.1)

zeta <- sqrt(CompI^2 + CompII^2)^2/4
zetareal <- sqrt(CompIreal^2 + CompIIreal^2)^2/4
zetaax <- sqrt((modeslin[1] + modeslin[5]*axval)^2 + (modeslin[6] + modeslin[10]*axval)^2)^2/4

Concentration <- sqrt((pi*zeta)/2)*exp(-zeta)*(besselI(zeta, 0) + besselI(zeta, 1))
Concentrationreal <- sqrt((pi*zetareal)/2)*exp(-zetareal)*(besselI(zetareal, 0) + besselI(zetareal, 1))
Concentrationax <- sqrt((pi*zetaax)/2)*exp(-zetaax)*(besselI(zetaax, 0) + besselI(zetaax, 1))
  
PredCirc <- atan2(CompII, CompI)
PredCircreal <- atan2(CompIIreal, CompIreal)


pdf("Pointing North Data/FigurelineSBSODreal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)

plot(Xreal, PredCircreal*(180/pi), ylim = c(-200,200), type = "l", xlab = "SBSOD",
     ylab = "Pointing Error", bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(-180, 0,180), at = c(seq(-180, 180, by = 60)),
     labels = c(seq(-180, 180, by = 60)), las = 2)
points(North$SBSOD, North$theta*(180/pi), cex=0.5)
points(axval,
       atan2(modeslin[6] + modeslin[10]*axval,
             modeslin[1] + modeslin[5]*axval)*(180/pi),
       pch = 0)

dev.off()


pdf("Pointing North Data/FigureConcentrationSBSODreal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)
plot(Xreal, Concentrationreal, type = "l", xlab = "SBSOD",ylab = "Concentration",
     bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(0, 0.5,1), at = c(seq(0, 1, by = 0.05)),
     labels = c(seq(0, 1, by = 0.05)), las = 2)
points(axval, Concentrationax, pch = 0)
dev.off()


#Age

X <- seq(-5, 5, 1)
Xreal <- seq(-1.68, 3.32, 0.1)
CompI  <- modeslin[1] + modeslin[2]*X
CompII <- modeslin[6] + modeslin[7]*X
CompIreal  <- modeslin[1] + modeslin[2]*Xreal
CompIIreal <- modeslin[6] + modeslin[7]*Xreal
axval <- hmode(Results$ax[1001:3000,1], 0.1)

zeta <- sqrt(CompI^2 + CompII^2)^2/4
zetareal <- sqrt(CompIreal^2 + CompIIreal^2)^2/4
zetaax <- sqrt((modeslin[1] + modeslin[2]*axval)^2 + (modeslin[6] + modeslin[7]*axval)^2)^2/4

Concentration <- sqrt((pi*zeta)/2)*exp(-zeta)*(besselI(zeta, 0) + besselI(zeta, 1))
Concentrationreal <- sqrt((pi*zetareal)/2)*exp(-zetareal)*(besselI(zetareal, 0) + besselI(zetareal, 1))
Concentrationax <- sqrt((pi*zetaax)/2)*exp(-zetaax)*(besselI(zetaax, 0) + besselI(zetaax, 1))

PredCirc <- atan2(CompII, CompI)
PredCircreal <- atan2(CompIIreal, CompIreal)


pdf("Pointing North Data/FigurelineAgereal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)

plot(Xreal, PredCircreal*(180/pi), ylim = c(-200,200), type = "l", xlab = "Age",
     ylab = "Pointing Error", bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(-180, 0,180), at = c(seq(-180, 180, by = 60)),
     labels = c(seq(-180, 180, by = 60)), las = 2)
points(North$SBSOD, North$theta*(180/pi), cex=0.5)
points(axval,
       atan2(modeslin[6] + modeslin[7]*axval,
             modeslin[1] + modeslin[2]*axval)*(180/pi),
       pch = 0)

dev.off()


pdf("Pointing North Data/FigureConcentrationAgereal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)
plot(Xreal, Concentrationreal, type = "l", xlab = "Age",ylab = "Concentration",
     bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(0, 0.5,1), at = c(seq(0, 1, by = 0.05)),
     labels = c(seq(0, 1, by = 0.05)), las = 2)
points(axval, Concentrationax, pch = 0)
dev.off()


#Experience

X <- seq(-5, 5, 1)
Xreal <- seq(-1.79, 2.21, 0.1)
CompI  <- modeslin[1] + modeslin[4]*X
CompII <- modeslin[6] + modeslin[9]*X
CompIreal  <- modeslin[1] + modeslin[4]*Xreal
CompIIreal <- modeslin[6] + modeslin[9]*Xreal
axval <- hmode(Results$ax[1001:3000,3], 0.1)

zeta <- sqrt(CompI^2 + CompII^2)^2/4
zetareal <- sqrt(CompIreal^2 + CompIIreal^2)^2/4
zetaax <- sqrt((modeslin[1] + modeslin[4]*axval)^2 + (modeslin[6] + modeslin[9]*axval)^2)^2/4

Concentration <- sqrt((pi*zeta)/2)*exp(-zeta)*(besselI(zeta, 0) + besselI(zeta, 1))
Concentrationreal <- sqrt((pi*zetareal)/2)*exp(-zetareal)*(besselI(zetareal, 0) + besselI(zetareal, 1))
Concentrationax <- sqrt((pi*zetaax)/2)*exp(-zetaax)*(besselI(zetaax, 0) + besselI(zetaax, 1))

PredCirc <- atan2(CompII, CompI)
PredCircreal <- atan2(CompIIreal, CompIreal)


pdf("Pointing North Data/FigurelineExpreal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)

plot(Xreal, PredCircreal*(180/pi), ylim = c(-200,200), type = "l", xlab = "Experience",
     ylab = "Pointing Error", bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(-180, 0,180), at = c(seq(-180, 180, by = 60)),
     labels = c(seq(-180, 180, by = 60)), las = 2)
points(North$SBSOD, North$theta*(180/pi), cex=0.5)
points(axval,
       atan2(modeslin[6] + modeslin[9]*axval,
             modeslin[1] + modeslin[4]*axval)*(180/pi),
       pch = 0)

dev.off()


pdf("Pointing North Data/FigureConcentrationExpreal.pdf", height = 5, width = 8,
    family = "CM Roman", pointsize = 12)
plot(Xreal, Concentrationreal, type = "l", xlab = "Experience",ylab = "Concentration",
     bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-3, 0, 3), at = c(seq(-3, 3, by = 1)),
     labels = c(seq(-3, 3, by = 1)), las = 2)
axis(2, xaxp = c(0, 0.5,1), at = c(seq(0, 1, by = 0.05)),
     labels = c(seq(0, 1, by = 0.05)), las = 2)
points(axval, Concentrationax, pch = 0)
dev.off()