source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
# loaded BC0.cls, BC0.dcor, BCBig.cls, BCBig.dcor

summary.plots(BC0.cls, BC0.dcor)
summary.plots(BCBig.cls, BCBig.dcor)


CLS <- BC0.cls
DCOR <- BC0.dcor
asGlyphs=FALSE; pad=FALSE; N=15; MIN=0.1; MAX=0.8; MOST=1; LEAST=0; DCOR.weight=2; useRaster=FALSE

R <- splom(BC0.cls, BC0.dcor)
R <- splom(BC0.cls, BC0.dcor, MAX=0.5, high.sat=FALSE)
R <- splom(BC0.cls, BC0.dcor, high.sat=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=FALSE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE, pad=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE, pad=TRUE, grid.col="#FF0000")

# EXAMPLE RASTER OUTPUT
scalar <- 10
width <- ncol(CLS)*scalar
height <- nrow(CLS)*scalar
png('test.raster.png', width=width, height=height, units="px", bg="white")
par(mar = rep(0, 4)) # set plot margins to 0 before drawing
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE, useRaster=TRUE, lwd=2)
dev.off()

# EXAMPLE VECTOR OUTPUT
width <- 7
height <- nrow(CLS)/ncol(CLS)*width
pdf('test.raster.pdf', width=width, height=height)
par(mar = rep(0, 4)) # set plot margins to 0 before drawing
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE)
dev.off()
