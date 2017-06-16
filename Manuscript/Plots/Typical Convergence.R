#Empty global environment
rm(list=ls())

#Load/source required packages/files
require(Rcpp)
sourceCpp('Paper 3 simulations/SimulationCode.cpp')
Res <- SimRegRCPP1P(50, 2, c(0,2), c(3,1), 0, 1, 5000, 1, 0, 500, 1000)


require(plotrix)
library(extrafont)
font_install("fontcm")
loadfonts(device = "win")
par(family = "fontcm")

png("Manuscript/Plots/TypicalHist.png", height = 6, width = 6,
    units = "in", family = "CM Roman", pointsize = 12, res = 1200)
par(mfrow = c(2,2))

hist(Res$B[1001:5000,1,9], breaks = 50, xlab = expression(beta[0]^I), main = "")
hist(Res$B[1001:5000,2,9], breaks = 50, xlab = expression(beta[1]^I), main = "")
hist(Res$B[1001:5000,3,9], breaks = 50, xlab = expression(beta[0]^II), main = "")
hist(Res$B[1001:5000,4,9], breaks = 50, xlab = expression(beta[1]^II), main = "")

dev.off()
