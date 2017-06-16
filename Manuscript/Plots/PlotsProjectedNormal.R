library(ggplot2)
library(dplyr)
library(grid)


integerTextCheck <- function(string) {
  toi <- strtoi(string)
  if (!is.na(toi)) return(toi) else return(string)
}


# Get a list of Projected Normal visualisations and tables.
getPNList <- function(b_0, b_1,
                      circResolution = 200,
                      xRange = c(-2, 2),
                      rangeExpand = 1.2,
                      nudgeExpand = 1,
                      getPlotRangeFromXRange = TRUE,
                      predCurve = FALSE,
                      predArrow = TRUE,
                      predDots = FALSE,
                      adda_x = TRUE,
                      adda_c = TRUE,
                      adda_c_dot = TRUE,

                      addXLabels = TRUE,
                      addXZeroLabel = TRUE,
                      addXEndpointLabels = TRUE,
                      xLabXis = TRUE,

                      addCircLines = TRUE,

                      circleColor = rgb(.8, .8, .8, .9),

                      inflectionShape = 1,
                      predLinePointShape = 1,
                      circlePointShape = 2,

                      inflectionSize    = 3,
                      predLinePointSize = 1,
                      circlePointSize   = 1,

                      xlim = NA,
                      ylim = NA,

                      BW = FALSE, # Display everything Black & White?
                      XthDegrees = TRUE, # Is the side plot in degrees
                      XthHalfCircle    = TRUE, # Display the full circle for xth
                      addCircPoints = TRUE, # Add points on the circle for x integers
                      labelBox = TRUE # Box around the x's

) {


  inflectionShape <- integerTextCheck(inflectionShape)
  predLinePointShape <- integerTextCheck(predLinePointShape)
  circlePointShape <- integerTextCheck(circlePointShape)

  outList <- list()

  # Circular variants of the predictors.
  if (b_1[1] != 0 & b_1[2] != 0) {
    a_x   <- -((b_0[1]*b_1[1] + b_0[2]*b_1[2]) / (b_1[1]^2 + b_1[2]^2))
    a_cI  <- b_0[1] + b_1[1] * a_x
    a_cII <- b_0[2] + b_1[2] * a_x
    a_c <- atan2(a_cII, a_cI) %% (2*pi)


  } else {
    a_x <- 0
    a_cI <- b_0[1]
    a_cII <- b_0[2]
    a_c <- atan2(a_cII, a_cI) %% (2*pi)
  }

  if (a_x != 0) {
    b_c <- tan(atan2(b_0[2], b_0[1])-a_c)/(-a_x)
  } else {
    b_c <- tan(atan2(b_0[2] + b_1[2], b_0[1] + b_1[1]) - a_c)/(1-a_x)
  }

  # THe shortest distance to the origin
  shortestOrigDist <- sqrt((b_0[1] + b_1[1]*a_x)^2 +
                             ((b_0[2] + b_1[2]*a_x)^2))


  if (shortestOrigDist < 2e-10) {
    b_c <- 1000000
    a_c <- a_c + pi/2

  }

  circLinesGeom  <- circPointGeom <- a_cPathGeom <- a_cLineGeom <-
    labelGeom <- a_cZeroAngleGeom <- geom_chosen_lty <- myArrow <-
    predDotGeom <- circLinesGeom <- circPointGeom <-
    a_xGeom <- a_xXthGeom <- a_c_dotGeom <- NULL


  ### PLOTTING

  # Some basic colors.
  if (!BW) {
    circLineColor  <- rgb(.3, .3, .4, .25)
    circPointColor <- rgb(.8, .3, .2, 1)
    predLineColor  <- rgb(.6, .6, .9, .8)
    predDotsColor  <- rgb(.7, .7, .9, 1)
    a_xDotColor    <- rgb(.9, .1, .6, 1)
    a_cDotColor    <- rgb(.8, .3, .6, 1)
  } else {
    circLineColor  <- rgb(.8, .8, .8, .9)#rgb(.3, .3, .3, .5)
    circPointColor <- rgb(.3, .3, .3, .5)#rgb(.8, .8, .8, .9)
    predLineColor  <- rgb(.6, .6, .6, .9)
    predDotsColor  <- rgb(.3, .3, .3, .5)#rgb(.7, .7, .7, .8)
    a_xDotColor    <- rgb(.1, .1, .1, .8)
    a_cDotColor    <- rgb(.6, .6, .6, .7)
  }


  # Limits for the base case.
  xmin <- min(-1, b_0[1]) * rangeExpand
  xmax <- max( 1, b_0[1]) * rangeExpand
  ymin <- min(-1, b_0[2]) * rangeExpand
  ymax <- max( 1, b_0[2]) * rangeExpand

  # Sequence of xvalues to test.
  xSeq <- seq(xRange[1], xRange[2], by = 1)

  # All the suggested points on the prediction line by xRange.
  fullPredLine <- data.frame(cos = b_0[1] + b_1[1]*xSeq,
                             sin = b_0[2] + b_1[2]*xSeq,
                             dir = atan2(b_0[2] + b_1[2]*xSeq,
                                         b_0[1] + b_1[1]*xSeq),
                             len = sqrt((b_0[1] + b_1[1]*xSeq)^2 +
                                          (b_0[2] + b_1[2]*xSeq)^2))
  row.names(fullPredLine) <- xSeq



  # Make limits larger to accomodate the whole xRange.
  if (getPlotRangeFromXRange) {
    xmin <- min(c(xmin, fullPredLine$cos)) * rangeExpand
    xmax <- max(c(xmax, fullPredLine$cos)) * rangeExpand
    ymin <- min(c(ymin, fullPredLine$sin)) * rangeExpand
    ymax <- max(c(ymax, fullPredLine$sin)) * rangeExpand

    # In this case we've made space for all values in xRange.
    predLine <- fullPredLine

  } else {

    # All the points that could be plotted with this plot size.
    predValuesInPlot <- (fullPredLine$cos > xmin & fullPredLine$cos < xmax &
                           fullPredLine$sin > ymin & fullPredLine$sin < ymax)

    # In this case, the predictor line might have to be made shorter.
    predLine <- fullPredLine[predValuesInPlot, ]
  }

  # Limits for the plot
  xl <- c(xmin, xmax)
  yl <- c(ymin, ymax)

  if (!is.na(xlim)) xl <- xlim
  if (!is.na(ylim)) yl <- ylim


  # Make a bare-bones base plot.
  theme_spartan <- theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  bp <- ggplot(data.frame()) + theme_spartan

  # The usual circular sequence, and its Euclidean counterpart.
  sq <- seq(0, 2 * pi, length.out = circResolution)
  xysq <- cbind(cos(sq), sin(sq))

  # Circle geom
  circ <- geom_path(aes(x = xysq[,1], y = xysq[,2]), col = circleColor)

  # Intercept point
  icpt <- geom_point(aes(x = b_0[1], y = b_0[2]))

  # Start building out
  outPlot <- bp + circ +
    xlim(xl) + ylim(yl) + coord_fixed() +
    xlab("") + ylab("")


  # Predictor line
  nPredPoints <- nrow(predLine)

  # Set the chosen parameters.
  if (predArrow) myArrow <- arrow(length = unit(0.03, "npc"), type="closed")
  if (predCurve) {geom_chosen_lty <- geom_curve} else {geom_chosen_lty <- geom_segment}

  predLineGeom <-  geom_chosen_lty(aes(x = predLine[-nPredPoints,1],
                                       y = predLine[-nPredPoints,2],
                                       xend = predLine[-1,1],
                                       yend = predLine[-1,2]),
                                   arrow = myArrow, color = predLineColor)

  predLineBGGeom <- geom_segment(aes(x = c()))

  if (predDots)  predDotGeom <- geom_point(aes(predLine$cos, predLine$sin),
                                           color = predDotsColor,
                                           size = predLinePointSize,
                                           shape = predLinePointShape)




  # Add labels for the x lines.
  if (addXLabels | addXZeroLabel | addXEndpointLabels) {


    # Set the desired x labels.
    if (addXLabels) {
      labeledPredLine <- predLine
    } else {

      predLineValsToBeLabeled <- rep(FALSE, nrow(predLine))

      if (addXZeroLabel & (0 %in% rownames(predLine))) {
        predLineValsToBeLabeled[rownames(predLine) == 0] <- TRUE
      }
      if (addXEndpointLabels) {
        predLineValsToBeLabeled[1] <- TRUE
        predLineValsToBeLabeled[nPredPoints] <- TRUE
      }

      labeledPredLine <- predLine[predLineValsToBeLabeled, ]
    }

    textAes <- aes(x = labeledPredLine$cos,
                   y = labeledPredLine$sin,
                   label = paste(ifelse(xLabXis, "x =", ""),
                                 row.names(labeledPredLine)))



    # Direction of the prediction line.
    predDirection  <- atan2(b_1[2], b_1[1]) %% (2*pi)

    # The nudge direction is ortogonal to the prediction direction.
    nudgeDirection <- (predDirection - pi/2) %% (2*pi)

    # Don't ask me why this works.
    xNudge <- sign(b_c) * cos(nudgeDirection)
    yNudge <- sign(b_c) * sin(nudgeDirection)

    if (labelBox) {
      labelGeom <- geom_label(mapping=textAes,
                              nudge_x = xNudge * nudgeExpand,
                              nudge_y = yNudge * nudgeExpand)
    } else {
      labelGeom <- geom_text(mapping=textAes,
                             size = 6,
                             nudge_x = xNudge * nudgeExpand,
                             nudge_y = yNudge * nudgeExpand)
    }
  }


  # Add lines to explain the effect on the circle.
  if (addCircLines | addCircPoints) {

    # Expand the lines for the points within the circle, to reach the end.
    circLineExpand <- 1/sapply(predLine$len, function(x) min(x, 1))

    circLineAes  <- aes(x = predLine$cos * circLineExpand,
                        y = predLine$sin * circLineExpand,
                        xend = 0, yend = 0)

    # Points on the circle.
    circPointAes <- aes(x = cos(predLine$dir)[predLine$len > 0],
                        y = sin(predLine$dir)[predLine$len > 0])

    if (addCircLines) circLinesGeom <- geom_segment(circLineAes,
                                                    color = circLineColor)
    if (addCircPoints) circPointGeom <- geom_point(circPointAes,
                                                   color = circPointColor,
                                                   size = circlePointSize,
                                                   shape = circlePointShape)
  }



  if (adda_x) a_xGeom <- geom_point(aes(a_cI, a_cII),
                                    color = a_xDotColor,
                                    size = inflectionSize,
                                    shape = inflectionShape)



  if (adda_c) {
    a_cCircSeq  <- seq(0, a_c, length.out = round(circResolution * a_c / (2*pi)))
    a_cXySeq    <-  cbind(cos(a_cCircSeq), sin(a_cCircSeq)) * 0.45

    a_cPathGeom <- geom_path(aes(x = a_cXySeq[,1], y = a_cXySeq[,2]),
                             color = a_cDotColor)

    a_cLineGeom <- geom_segment(aes(x = cos(a_c), y = sin(a_c),
                                    xend = 0, yend = 0), color = a_cDotColor)
    a_cZeroAngleGeom <- geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0),
                                     color = rgb(0, 0, 0, .2))
  }

  if (adda_c_dot) {
    a_c_dotGeom <- geom_point(aes(cos(a_c), sin(a_c)),
                              color = a_xDotColor,
                              size = inflectionSize,
                              shape = inflectionShape)
  }




  predPlot <- outPlot + a_cZeroAngleGeom + predLineGeom + predDotGeom +
    circLinesGeom + circPointGeom + labelGeom +
    a_cPathGeom + a_cLineGeom + a_xGeom + a_c_dotGeom




  # A table of predicted values.
  predTable <- data.frame(x = row.names(predLine),
                          "CircPrediction" = predLine$dir %% (2*pi))

  # A table of parameter values.
  paramTable <- t(data.frame(b_0_I = b_0[1], b_0_II = b_0[2],
                             b_1_I = b_1[1], b_1_II = b_1[2],
                             a_x = a_x, a_c = a_c, b_c = b_c))
  colnames(paramTable) <- "Value"


  predFun <- function(x) ((a_c + atan(b_c * (x - a_x)) - a_c + pi) %% (2*pi)) + a_c - pi


  # If we have to display in degrees, use the correct transformation factor.
  if (XthDegrees) {
    degFac <- 180/pi
  } else {
    degFac <- 1
  }

  if (XthHalfCircle) {
    XthYlim <- ylim((a_c - pi/2)*degFac, (a_c + pi/2)*degFac)
  } else {
    XthYlim <- NULL
  }


  if(adda_x) a_xXthGeom <- geom_point(aes(x = a_x, y = a_c * degFac),
                                      color = a_xDotColor,
                                      size = inflectionSize,
                                      shape = inflectionShape)




  xthPlot <- ggplot(data.frame()) +
    xlim(c(xRange[1] - 0.9, xRange[2] + 0.9)) + XthYlim +
    stat_function(aes(x = xRange), fun = function(x) predFun(x) * degFac) +
    geom_point(aes(x = xSeq, y = predFun(xSeq) * degFac),
               color = circPointColor,
               shape = circlePointShape) +
    a_xXthGeom +
    xlab("Predictor value x") + ylab("Circular predicted value") +
    theme_minimal()



  list(predPlot = predPlot)

}

#####Create the Plot


Plot2 <- getPNList(c(0.5, -0.5), c(0.15, 0.27),
                   circResolution = 200,
                   xRange = c(-3, 3),
                   rangeExpand = 1.1,
                   nudgeExpand = -0.3,
                   getPlotRangeFromXRange = TRUE,
                   predCurve = FALSE,
                   predArrow = FALSE,
                   predDots = TRUE,
                   adda_x = TRUE,
                   adda_c = FALSE,
                   adda_c_dot = TRUE,

                   addXLabels = FALSE,
                   addXZeroLabel = TRUE,
                   addXEndpointLabels = TRUE,
                   xLabXis = FALSE,

                   addCircLines = FALSE,

                   inflectionShape = 21,
                   predLinePointShape = 1,
                   circlePointShape = 16,

                   inflectionSize    = 0,
                   predLinePointSize = 2,
                   circlePointSize   = 2,

                   xlim = c(-1.4, 1.4),
                   ylim = c(-1.4, 1.4),

                   BW = TRUE, # Display everything Black & White?
                   XthDegrees = TRUE, # Is the side plot in degrees
                   XthHalfCircle    = TRUE, # Display the full circle for xth
                   addCircPoints = TRUE, # Add points on the circle for x integers
                   labelBox = FALSE # Box around the x's

)

Plot3 <- getPNList(c(0.5, -0.5), c(0.15, -0.15),
                   circResolution = 200,
                   xRange = c(-3, 3),
                   rangeExpand = 1.1,
                   nudgeExpand = -0.3,
                   getPlotRangeFromXRange = TRUE,
                   predCurve = FALSE,
                   predArrow = FALSE,
                   predDots = TRUE,
                   adda_x = TRUE,
                   adda_c = FALSE,
                   adda_c_dot = TRUE,

                   addXLabels = FALSE,
                   addXZeroLabel = TRUE,
                   addXEndpointLabels = TRUE,
                   xLabXis = FALSE,


                   addCircLines = FALSE,

                   inflectionShape = 21,
                   predLinePointShape = 1,
                   circlePointShape = 16,

                   inflectionSize    = -0.1,
                   predLinePointSize = 2,
                   circlePointSize   = 2,

                   xlim = c(-1.4, 1.4),
                   ylim = c(-1.4, 1.4),

                   BW = TRUE, # Display everything Black & White?
                   XthDegrees = TRUE, # Is the side plot in degrees
                   XthHalfCircle    = TRUE, # Display the full circle for xth
                   addCircPoints = TRUE, # Add points on the circle for x integers
                   labelBox = FALSE # Box around the x's

)


require(plotrix)
library(extrafont)
font_install("fontcm")
loadfonts()
loadfonts(device = "win")
par(family = "LM Roman 10")

png("Manuscript/Plots/Location.png", width = 5, height = 5, family = "LM Roman 10", pointsize = 20,
    units = "in", res = 1200)
Plot2
dev.off()


png("Manuscript/Plots/Accuracy.png", width = 5, height = 5, family = "LM Roman 10", pointsize = 20,
    units = "in", res = 1200)
Plot3
dev.off()



