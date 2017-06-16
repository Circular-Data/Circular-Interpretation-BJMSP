require(plotrix)
require(shape)
library(extrafont)
font_install("fontcm")
loadfonts()
loadfonts(device = "win")
par(family = "LM Roman 10")


png("Manuscript/Plots/SSDO.png", width=7, height=7, units = "in", res = 1200, family = "LM Roman 10")

x <- seq(-14, 14, 1)
y <- seq(-14, 14, 1)

cutoff <- atan2(3,3)
ac1 <- atan2(-3, 3)
ac2 <- atan2(3, -3)

plot(1, type="n", xlab="", ylab="", yaxt = "n", xaxt= "n", xlim=c(-10, 10), ylim=c(-10, 10), frame.plot = FALSE)
#axis(side=2)
#axis(side=1)
polygon(c(-10, -10, 10), c(-10, 10, 10), col = "grey", border = NA)
draw.circle(0, 0, 3, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1)

#draw lines
points(x, y, type = "l", lty = 2)
Arrowhead(5,5,45)
#intersections
points(cos(cutoff)*3+0.1, sin(cutoff)*3+0.1, pch = 0)
points(cos(cutoff2)*3-0.1, sin(cutoff2)*3-0.1, pch = 0)
points(x, y+5, type = "l")
Arrowhead(3,8,45)
#ac1
points(cos(ac2)*3-0.1, sin(ac2)*3+0.1, pch = 16)
text(cos(ac2)*3+0.3, sin(ac2)*3-0.3, label = expression(a[c]))
points(x, y-5, type = "l")
Arrowhead(7,2,45)
#ac2
points(cos(ac1)*3+0.1, sin(ac1)*3-0.1, pch = 16)
text(cos(ac1)*3-0.3, sin(ac1)*3+0.3, label = expression(a[c]))

text(3.4, 0, label = 0)
text(-3.4, 0.4, label = expression(pi))
text(-3.6, -0.4, label = expression(-pi))
text(-9,-7, label = expression(paste("(",beta[1]^I,"x",",",~beta[1]^II,"x",")")))
text(-7,3, label = expression(paste( "(",beta[0]^I + beta[1]^I, "x",",", ~beta[0]^II + beta[1]^II,"x",")")))
text(3,-7, label = expression(paste("(",-beta[0]^I + beta[1]^I, "x",",", ~-beta[0]^II + beta[1]^II,"x",")")))

#Close the device
dev.off()





png("Manuscript/Plots/SSDO2.png", width=7, height=7, units = "in", res = 1200, family = "LM Roman 10")

x <- seq(-14, 14, 1)
y <- seq(-14, 14, 1)

cutoff <- atan2(3,3)
ac2 <- atan2(-3, 3)
ac <- atan2(3, -3)

plot(1, type="n", xlab="", ylab="", yaxt = "n", xaxt= "n", xlim=c(-10, 10), ylim=c(-10, 10), frame.plot = FALSE)
#axis(side=2)
#axis(side=1)
polygon(c(-10,-10, 10, 10), c(0, 10, 10, 0), col = "grey", border = NA)
draw.circle(0, 0, 3, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1)

#draw lines
points(x, rep(0,29) , type = "l", lty = 2)
points(x, (rep(5,29)), type = "l")
Arrowhead(5,0,0)
Arrowhead(5,5,0)
#intersections
points(cos(ac-cutoff)*3, sin(ac-cutoff)*3+0.2, pch = 2)

text(3.4, 0, label = 0)
text(-3.4, 0.4, label = expression(pi))
text(-3.6, -0.4, label = expression(-pi))

#Close the device
dev.off()





png("Manuscript/Plots/SSDO3.png", width=7, height=7, units = "in", res = 1200, family = "LM Roman 10")

x <- seq(-14, 14, 1)
y <- seq(-14, 14, 1)

cutoff <- atan2(3,3)
ac2 <- atan2(-3, 3)
ac <- atan2(3, -3)

plot(1, type="n", xlab="", ylab="", yaxt = "n", xaxt= "n", xlim=c(-10, 10), ylim=c(-10, 10), frame.plot = FALSE)
#axis(side=2)
#axis(side=1)
polygon(c(-10,-10, 10, 10), c(0, 10, 10, 0), col = "grey", border = NA)
draw.circle(0, 0, 3, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1)

#draw lines
points(x, rep(0,29) , type = "l", lty = 2)
points(x, (rep(-5,29)), type = "l")
Arrowhead(5,0,0)
Arrowhead(5,-5,0)
#intersections
points(cos(ac2-cutoff)*3, sin(ac2-cutoff)*3-0.3, pch = 2)

text(3.4, 0, label = 0)
text(-3.4, 0.4, label = expression(pi))
text(-3.6, -0.4, label = expression(-pi))


#Close the device
dev.off()