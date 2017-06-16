library(extrafont)
font_import() #only run once
loadfonts(device = "win")
font_install("LM Roman 10")
font_install("fontcm")
par(family = "LM Roman 10")
set.seed(101)

Xreal <- rnorm(100, -2, 1.5)

CompI  <- 1 + 0.3*Xreal + rnorm(100, 0.1,1)
CompII <- 2 + -0.4*Xreal + rnorm(100, 0.1,1)

Outreal <- atan2(CompII, CompI)

X <- seq(-5, 5, 0.1)

predCompI  <- 1.1 + 0.25*X
predCompII <- 2.1 + -0.4*X

PredCirc <- atan2(predCompII, predCompI)

axval <- ((1.1*0.25)+(2.1*-0.4))/(1.1^2 + 2.1^2)
acval <- atan2(2.1 + -0.4*axval, 1.1 + 0.25*axval)

png("Manuscript/Plots/CircRegLine.png", height = 5, width = 8,
    family = "LM Roman 10", units = "in", pointsize = 12, res = 1200)

plot(X, PredCirc*(180/pi), ylim = c(-200,200), type = "l", xlab = "Predictor",
     ylab = "Circular Outcome", bty = 'n', yaxt = "n", xaxt = "n")
axis(1, xaxp = c(-5, 0, 5), at = c(seq(-5, 5, by = 1)),
     labels = c(seq(-5, 5, by = 1)), las = 2)
axis(2, xaxp = c(-180, 0,180), at = c(seq(-180, 180, by = 60)),
     labels = c(seq(-180, 180, by = 60)), las = 2)
points(Xreal, Outreal*(180/pi), cex=0.5)
points(axval,
       acval*(180/pi),
       pch = 0)

dev.off()
