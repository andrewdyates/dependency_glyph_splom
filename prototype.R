##
## * Plot color glyphs enumeration
## * Plot color dcor
## Plot color glyph matrix
## Plot scaled colored glyph matrix
## return top dCOR per class
##
library("RColorBrewer")
library(modeest)
library("gplots")

BC0.cls <- as.matrix(read.table("data/bc0.cls.tab", header=TRUE, sep=" ", row.names=1))
BC0.dcor <- as.matrix(read.table("data/bc0.dcor.tab", header=TRUE, sep=" ", row.names=1))
BCBig.cls <- as.matrix(read.table("data/bcBig.cls.tab.gz", header=TRUE, sep=" ", row.names=1))
BCBig.dcor <- as.matrix(read.table("data/bcbig.dcor.tab.gz", header=TRUE, sep=" ", row.names=1))
GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")
GRID.COL <-  "#cdd7e6"
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
BC0.cls <- renumerate(BC0.cls)
BCBig.cls <- renumerate(BCBig.cls)
save(BC0.cls, BC0.dcor, BCBig.cls, BCBig.dcor, file="BC.RData")
# load("BC.RData")

CLS <- as.matrix(BC0.cls)
DCOR <- as.matrix(BC0.dcor)
# ==============================
# COMBINED GLYPH AND DCOR SPLOM
# --------------------
# Boolean Class distance
D.cls.r <- dist(CLS)
D.cls.c <- dist(t(CLS))
# dCOR distance. WARNING: cor on columns, dist on rows!
D.dcor.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)
D.dcor.c <- as.dist(1-cor(DCOR, method="pearson")) + dist(t(DCOR))
# Row and column ordering based on mean dCOR
Rowv <- rowMeans(DCOR, na.rm = TRUE)
Colv <- colMeans(DCOR, na.rm = TRUE)

# Hierarchical cluster on both distance correlation and class distance
Rhclust <- as.dendrogram(hclust(D.dcor.r*3+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
Chclust <- as.dendrogram(hclust(D.dcor.c*3+D.cls.c, method="average"))
Chclust <- reorder(Chclust, Colv)

heatmap_cols <- c("#ffffff", "#a00d42", "#d7424c", "#eb6532", "#000000", "#40a185", "#2688bf", "#5b51a5")
heatmap.2(expand.cls(as.matrix(CLS)),
  col=heatmap_cols, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=0:8-0.5,
  key=TRUE, symkey=FALSE, trace="none",
  Rowv=Rhclust, Colv=Chclust
);
# ------------------------------



CLS <- as.matrix(BC0.cls)
DCOR <- as.matrix(BC0.dcor)
# ==============================
# DCOR
# ------------------------------
MIN.DCOR<-0.08
MAX.DCOR<-0.8
heatmap_breaks_dcor <- seq(MIN.DCOR,MAX.DCOR,0.01)
heatmap_cols_dcor <- rev(colorRampPalette(brewer.pal(8,"RdYlBu"))(length(heatmap_breaks_dcor)-1))
# Hierarchical cluster on both distance correlation and class distance
# dCOR distance. WARNING: cor on columns, dist on rows!
D.dcor.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)
D.dcor.c <- as.dist(1-cor(DCOR, method="pearson")) + dist(t(DCOR))
# Row and column ordering based on mean dCOR
Rowv <- rowMeans(DCOR, na.rm = TRUE)
Colv <- colMeans(DCOR, na.rm = TRUE)
Rhclust <- as.dendrogram(hclust(D.dcor.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
Chclust <- as.dendrogram(hclust(D.dcor.c, method="average"))
Chclust <- reorder(Chclust, Colv)

heatmap.2(as.matrix(DCOR), 
  col=heatmap_cols_dcor, ylab="CpG", xlab="mRNA", symm=FALSE, breaks=heatmap_breaks_dcor,
  key=TRUE, symkey=FALSE, trace="none",
  Rowv=Rhclust, Colv=Chclust
);
# ------------------------------



CLS <- as.matrix(BC0.cls)
DCOR <- as.matrix(BC0.dcor)
# ==============================
# COMBINED GLYPH SPLOM, NO SCALING
# --------------------

### Expand CLS matrix into 2x2 glyphs.
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

G <- expand.cls(CLS)
Gb <- expand.cls(CLS, pad=TRUE)

## Plot enumerated glyph matrix as image. Works on both expanded and compressed glyph matrices.
## Draw grid parameter:
##   0: do not draw grid
##   1: draw grid every square
##   2: draw grid every other square (for 2x2 glyphs)
draw.glyphs <- function(G, col=GLYPH.COLS, draw.grid=0, grid.col=GRID.COL) {
  ## "Lower edge" of color bins. Bin is [i,i+1) of `breaks`, indexed from 1.
  ##    "0" goes into the first bin: [-0.5, 0.5)
  ##    "1" goes into second bin:    [0.5, 1.5)
  ##    ...
  ##    "7" goes into ninth bin (eight bins, plus lowest bound)  [7.5, 8.5)
  breaks <- 0:length(GLYPH.COLS)-0.5
  ## Convert matrix into "image" so that top left corners are aligned
  ##   `image` plots transpose of matrix
  ##   y-axis ploted from low to high index (from bottom left rather than top left corner)
  Img <- t(G)[,seq(nrow(G),1,-1)]
  w<-ncol(G); h<-nrow(G)
  image(1:w, 1:h, Img, col=col, breaks=breaks,
    axes=FALSE, xlab="", ylab="")
  if (draw.grid == 1) {
    for(i in 1:w)
      abline(v=i-0.5, untf=FALSE, col=grid.col)
    for(j in 1:h)
      abline(h=j-0.5, untf=FALSE, col=grid.col)
  } else if (draw.grid == 2) {
    for(i in 1:(w/2))
      abline(v=i*2-0.5, untf=FALSE, col=grid.col)
    for(j in 1:(h/2))
      abline(h=j*2-0.5, untf=FALSE, col=grid.col)
  }
}

## ORDER THE MATRIX
## --------------------
# Boolean Class distance
CLS <- BC0.cls
D.cls.r <- dist(CLS)
D.cls.c <- dist(t(CLS))
Rowv <- rowMeans(CLS, na.rm = TRUE)
Colv <- colMeans(CLS, na.rm = TRUE)

Rhclust <- as.dendrogram(hclust(D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
Chclust <- as.dendrogram(hclust(D.cls.c, method="average"))
Chclust <- reorder(Chclust, Colv)

# reorder rows and columns based on clustering
rowInd <- order.dendrogram(Rhclust)
colInd <- order.dendrogram(Chclust)

CLS <- CLS[rowInd, colInd]
draw.glyphs(CLS)
draw.glyphs(expand.cls(CLS))
