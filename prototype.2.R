## Scale glyph SPLOM intensity by dCOR
##
source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
# loaded BC0.cls, BC0.dcor, BCBig.cls, BCBig.dcor

# colorRampPalette

CLS  <- BC0.cls
DCOR <- BC0.dcor
R <- get.cls.order(CLS)
rowInd <- order.dendrogram(R$Rhclust)
colInd <- order.dendrogram(R$Chclust)
CLS.ord <- CLS[rowInd, colInd]
DCOR.ord <- DCOR[rowInd, colInd]

# OK, can we scale the intensities?
# make n=10 levels per color, interpolate from white


GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")


# Row is level (1 is white, N is color); Col is glyph in GLYPH.COLS order


N <- 10
## First interpolation level is white for values below MIN.
## Thus there are N-1 tints, plus white, for each range.
## We assign values of at least MAX the untinted color.
COLOR.M <- sapply(GLYPH.COLS, function(color) colorRampPalette(c("#ffffff", color))(N))

MIN <- 0.08
MAX <- 0.8

# Mark class ordinals using a decimal (less than 1). Order by increasing dep strength.

## "Lower edge" of color bins. Bin is [i,i+1) of `breaks`, indexed from 1.
##    "0" goes into the first bin: [-0.5, 0.5)
##    "1" goes into second bin:    [0.5, 1.5)
##    ...
##    "7" goes into ninth bin (eight bins, plus lowest bound)  [7.5, 8.5)

th <- (MAX-MIN)/(N-2)   ## bin for above and below bin threshold
breaks <- sapply(0:(N-2), function(i) MIN+i*th)
breaks <- c(0, breaks, 1)

## Select greatest break less than x
f <- function(x) min(tail(breaks[breaks<=0.7],1), MAX)
## Generate offset fraction from DCOR histogram.
OFFSET <- apply(DCOR, c(1,2), f)


# CLS+OFFSET ...
# convert COLOR.M to vector
# draw.glyphs(CLS.ord)
## ------------------------------
# generate color breaks from 0:7
#
# how to get concatenated matrix?
# HOW TO RESHAPE MATRIX TO A VECTOR?
B <- sapply(0:7, function(x) rep(x,N)+breaks[1:N])

color.breaks 
w<-ncol(G); h<-nrow(G)
image(1:w, 1:h, Img, col=col, breaks=breaks,
  axes=FALSE, xlab="", ylab="", useRaster=useRaster)

Img <- t(G)[,seq(nrow(G),1,-1)]



