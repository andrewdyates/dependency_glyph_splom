library("RColorBrewer")
library("gplots")

GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")
GRID.COL <-  "#ffffff"
# ---------------
# Re-enumerate glyph symbols to put ND in the center, 0 as bg
# Input enumeration:
##   'UNL': 0,
##   'HIH': 1,
##   'PC': 2,
##   'LIL': 3,
##   'HIL': 4,
##   'NC': 5,
##   'LIH': 6,
##   'NA': 7,
# ------------------------------
# Output enumeration:
##   'NA': 0,      white   => #ffffff
##   'HIH': 1,     teal    => #40a185
##   'PC': 2,      blue    => #2688bf
##   'LIL': 3,     purple  => #5b51a5
##   'UNL': 4,     black   => #000000
##   'HIL': 5,     magenta => #a00d42
##   'NC': 6,      red     => #d7424c
##   'LIH': 7,     orange  => #eb6532
##   [pad]: 8,     white   => #ffffff
renumerate <- function(BC) {
  BC[BC==0] <- -1
  BC[BC==7] <- 0
  BC[BC>=4 & BC<=6] <- BC[BC>=4 & BC<=6]+1
  BC[BC==-1] <- 4
  BC
}


## Extend heatmap.2 to generate a heatmap with decent colors and clustering.
## --------------------
## EXAMPLE:
## H <- heatmap.3(BC0.Dcor, ylab="CpG", xlab="mRNA", symm=FALSE)
heatmap.3 <- function(M, MIN=0.08, MAX=0.8, cols=brewer.pal(8,"RdYlBu"), ...) {
  heatmap_breaks <- seq(MIN,MAX,0.01)
  heatmap_cols <- rev(colorRampPalette(cols)(length(heatmap_breaks)-1))
  # distance as sum of euclidean and correlation. WARNING: cor on columns, dist on rows!
  D.r <- as.dist(1-cor(t(M), method="pearson")) + dist(M)
  D.c <- as.dist(1-cor(M, method="pearson")) + dist(t(M))
  # Row and column ordering based on mean values
  Rowv <- rowMeans(M, na.rm = TRUE)
  Colv <- colMeans(M, na.rm = TRUE)
  Rhclust <- as.dendrogram(hclust(D.r, method="average"))
  Rhclust <- reorder(Rhclust, Rowv)
  Chclust <- as.dendrogram(hclust(D.c, method="average"))
  Chclust <- reorder(Chclust, Colv)

  heatmap.2(as.matrix(M),
    col=heatmap_cols, breaks=heatmap_breaks,
    key=TRUE, symkey=FALSE, trace="none",
    Rowv=Rhclust, Colv=Chclust, ...
  )
}


## Cluster an enumerated boolean dependency class matrix.
## --------------------
# Helper function to return class ordering as Row / Column dendrogram objs.
get.cls.order <- function(CLS) {
  R = list()
  D.cls.r <- dist(CLS)
  D.cls.c <- dist(t(CLS))
  Rowv <- rowMeans(CLS, na.rm = TRUE)
  Colv <- colMeans(CLS, na.rm = TRUE)
  R$Rhclust <- as.dendrogram(hclust(D.cls.r, method="average"))
  R$Rhclust <- reorder(R$Rhclust, Rowv)
  R$Chclust <- as.dendrogram(hclust(D.cls.c, method="average"))
  R$Chclust <- reorder(R$Chclust, Colv)
  R
}
# Reorder class matrix by clustering.
order.cls <- function(CLS) {
  R <- get.cls.order(CLS)
  rowInd <- order.dendrogram(R$Rhclust)
  colInd <- order.dendrogram(R$Chclust)
  CLS[rowInd, colInd]
}


### Expand CLS matrix into 2x2 glyphs.
## --------------------
# 0,  1,2,3,  4,  5,6,7
to.glyph <- function(c) {
  r <- NaN
  if(c == 0)  # NA    (no significant dependency)
    r <- matrix(c(0,0,0,0), nrow=2)
  if(c == 1)  # HIH   (high x implies high y)
    r <- matrix(c(1,1,1,0), nrow=2)
  if(c == 2)  # PC    (positive correlation)
    r <- matrix(c(0,2,2,0), nrow=2)
  if(c == 3)  # LIL   (low x implies low y)
    r <- matrix(c(0,3,3,3), nrow=2)
  if(c == 4)  # UNL   (unspecified non-linear)
    r <- matrix(c(4,4,4,4), nrow=2)
  if(c == 5)  # HIL   (high x implies low y)
    r <- matrix(c(5,5,0,5), nrow=2)
  if(c == 6)  # NC    (negative correlation)   
    r <- matrix(c(6,0,0,6), nrow=2)
  if(c == 7)  # LIH   (low x implies low y)   
    r <- matrix(c(7,0,7,7), nrow=2)
  r
}


## Given a class matrix, construct a glyph matrix.
## --------------------
expand.cls <- function(CLS, pad=FALSE) {
  if(!pad) {
    G <- mat.or.vec(nrow(CLS)*2, ncol(CLS)*2)
    for(i in 0:(nrow(CLS)-1))
      for(j in 0:(ncol(CLS)-1))
        G[(i*2+1):(i*2+2),(j*2+1):(j*2+2)] <- to.glyph(CLS[i+1,j+1])
  } else {
    G <- mat.or.vec(nrow(CLS)*3+1, ncol(CLS)*3+1)
    G[,] <- 8 # fill with background enumeration
    for(i in 0:(nrow(CLS)-1))
      for(j in 0:(ncol(CLS)-1))
        G[(i*3+2):(i*3+3),(j*3+2):(j*3+3)] <- to.glyph(CLS[i+1,j+1])
  }
  G
}


## Plot enumerated glyph matrix as image. Works on both expanded and compressed glyph matrices.
## --------------------
## Draw grid parameter:
##   0: do not draw grid
##   1: draw grid every square
##   2: draw grid every other square (for 2x2 glyphs)
draw.glyphs <- function(G, col=GLYPH.COLS, grid=0, grid.col=GRID.COL, useRaster=FALSE, lwd=1) {
  ## "Lower edge" of color bins. Bin is (i,i+1] of `breaks`, indexed from 1.
  ##    "0" goes into the first bin: (-0.5, 0.5]
  ##    "1" goes into second bin:    (0.5, 1.5]
  ##    ...
  ##    "7" goes into ninth bin (eight bins, plus lowest bound)  (7.5, 8.5]
  breaks <- 0:length(GLYPH.COLS)-0.5
  ## Convert matrix into "image" so that top left corners are aligned
  ##   `image` plots transpose of matrix
  ##   y-axis ploted from low to high index (from bottom left rather than top left corner)
  Img <- t(G)[,seq(nrow(G),1,-1)]
  w<-ncol(G); h<-nrow(G)
  image(1:w, 1:h, Img, col=col, breaks=breaks,
    axes=FALSE, xlab="", ylab="", useRaster=useRaster)
  ## Add vector grid.
  if (grid != 0) {
    os <- grid-0.5
    if(grid==3) os<-os-0.5
    for(i in 0:(w/grid)+1)
      abline(v=i*grid-os, untf=FALSE, col=grid.col, lwd=lwd)
    for(j in 0:(h/grid)+1)
      abline(h=j*grid-os, untf=FALSE, col=grid.col, lwd=lwd)
  }
}


## Pixel-perfect plotting as png.
## --------------------
## EXAMPLE: plot.pix(expand.cls(CLS.ord), grid=2, scale=4, grid.col="white")
plot.pix <- function(G, fname="cls.px.plot.png", scale=1, ...) {
  fname <- paste0(sub("\\.png$", "", fname),".png")
  lwd=max(log(scale),1)
  width <- ncol(G)*scale
  height <- nrow(G)*scale
  png(fname, width=width, height=height, units="px", bg="white")
  par(mar = rep(0, 4)) # set plot margins to 0 before drawing
  draw.glyphs(G, useRaster=TRUE, lwd=lwd, ...)
  dev.off()
  # Return size of image saved in pixels.
  c(width, height)
}


## Scaled vector plotting as pdf. If height is NULL, scale proportionately with width.
## --------------------
## EXAMPLE: plot.vtr(expand.cls(CLS.ord), grid=2, grid.col="white")
plot.vtr <- function(G, width=7, height=NULL, fname="cls.plot.pdf", ...) {
  fname <- paste0(sub("\\.pdf$", "", fname),".pdf")
  scale <- width / ncol(G)
  if(is.null(height)) {
    height <- nrow(G)*scale
  }
  lwd=max(log(scale*72),1)
  pdf(fname, width=width, height=height)
  par(mar = rep(0, 4)) # set plot margins to 0 before drawing
  draw.glyphs(G, ...)
  dev.off()
  # Return size of image saved in inches.
  c(width, height)
}
