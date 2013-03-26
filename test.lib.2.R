source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
# loaded BC0.cls, BC0.dcor, BCBig.cls, BCBig.dcor
pdf("gibbs.BC.summary.plots.BC0.pdf")
summary.plots(BC0.cls, BC0.dcor)
dev.off()
pdf("gibbs.BC.summary.plots.BCBig.pdf")
summary.plots(BCBig.cls, BCBig.dcor)
dev.off()

# Figures for IEEE; BC0 (91,62)
pdf("gibbs.BC0.gsplom.pdf", width=10, height=14.68)
R <- splom(BC0.cls, BC0.dcor, MIN=0.1)
splom(BC0.cls, BC0.dcor, MIN=0.1, asGlyphs=T)
splom(BC0.cls, BC0.dcor, MIN=0.1, asGlyphs=T, lwd=10)
dev.off()
pdf("gibbs.BC0.gsplom.dendros.pdf", width=20, height=8)
plot(R$Rhclust, main="Row Dendrogram")
plot(R$Chclust, main="Column Dendrogram")
dev.off()

scalar <- 1
width <- ncol(BC0.cls)*scalar
height <- nrow(BC0.cls)*scalar
png('gibbs.BC0.gsplom.raster.1px.png', width=width, height=height, units="px", bg="white")
par(mar = rep(0, 4)) 
R <- splom(BC0.cls, BC0.dcor, asGlyphs=F, useRaster=TRUE, lwd=0, MIN=0.1, draw.labs=F)
dev.off()

png('gibbs.BC0.gsplom.raster.2px.png', width=width*2, height=height*2, units="px", bg="white")
par(mar = rep(0, 4)) 
R <- splom(BC0.cls, BC0.dcor, asGlyphs=T, useRaster=TRUE, lwd=0, MIN=0.1, draw.labs=F)
dev.off()
# ----------------------------------------
# Figures for IEEE; BCBig (1004,304)
pdf("gibbs.BCBig.gsplom.pdf", width=10, height=30)
R <- splom(BCBig.cls, BCBig.dcor, MIN=0.1)
qr <- order.dendrogram(R$Rhclust)
qc <- order.dendrogram(R$Chclust)
splom(BCBig.cls[qr,qc], BCBig.dcor[qr,qc], MIN=0.1, asGlyphs=T, reorder=F)
dev.off()
pdf("gibbs.BCBig.gsplom.dendros.pdf", width=20, height=8)
plot(R$Rhclust, main="Row Dendrogram")
plot(R$Chclust, main="Column Dendrogram")
dev.off()

scalar <- 1
width <- ncol(BC0.cls)*scalar
height <- nrow(BC0.cls)*scalar
png('gibbs.BCBig.gsplom.raster.1px.png', width=width, height=height, units="px", bg="white")
par(mar = rep(0, 4)) 
R <- splom(BCBig.cls, BCBig.dcor, asGlyphs=F, useRaster=TRUE, lwd=0, MIN=0.1, draw.labs=F)
dev.off()

png('gibbs.BCBig.gsplom.raster.2px.png', width=width*2, height=height*2, units="px", bg="white")
par(mar = rep(0, 4)) 
R <- splom(BCBig.cls, BCBig.dcor, asGlyphs=T, useRaster=TRUE, lwd=0, MIN=0.1, draw.labs=F)
dev.off()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CLS <- BC0.cls
DCOR <- BC0.dcor


R <- splom(BC0.cls, BC0.dcor)
R <- splom(BC0.cls, BC0.dcor, MAX=0.5, high.sat=FALSE)
R <- splom(BC0.cls, BC0.dcor, high.sat=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=FALSE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE, pad=TRUE)
R <- splom(BC0.cls, BC0.dcor, asGlyphs=TRUE, pad=TRUE, grid.col="#FF0000")

R <- splom(BC0.cls)
R <- splom(BC0.cls, asGlyphs=TRUE)
R <- splom(BC0.cls, asGlyphs=TRUE, grid.col="#FF0000")
R <- splom(BC0.cls, asGlyphs=TRUE, pad=TRUE, grid.col="#FF0000")

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
