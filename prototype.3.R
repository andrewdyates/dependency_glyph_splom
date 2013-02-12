## Generate glyphs without coloring
## Custom Color Scales
source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
load("data/methyl_mrna.RData")
CLS <- BC0.cls
DCOR <- BC0.dcor

ROW.VARS = methyl[which(rownames(methyl) %in% rownames(CLS)),]
COL.VARS = mRNA[which(rownames(mRNA) %in% colnames(CLS)),]

# Get pearson's correlation
V <- cor(t(ROW.VARS), t(COL.VARS))
heatmap.3(V, MIN=-1, MAX=1)
heatmap.3(sqrt(DCOR), MIN=0, MAX=1)
DELTA <- sqrt(DCOR)-abs(V)
heatmap.3(DELTA, MIN=min(DELTA), MAX=max(DELTA))

R <- get.order.cls.dcor(CLS, DCOR)
rowInd <- order.dendrogram(R$Rhclust)
colInd <- order.dendrogram(R$Chclust)
DELTA.ord <- DELTA[rowInd, colInd]
CLS.ord   <- CLS[rowInd, colInd]
DCOR.ord   <- DCOR[rowInd, colInd]

# Delta difference by glyph ordering
heatmap.3(DELTA.ord, MIN=min(DELTA), MAX=max(DELTA), preserve=TRUE)

# Expand matrix into glyphs; display

to.glyph <- function(z, bg=NA) {
  r <- NaN
  if(z >= 0 && z < 1)  # NA    (no significant dependency)
    r <- matrix(c(bg,bg,bg,bg), nrow=2)
  if(z >= 1 && z < 2)  # HIH   (high x implies high y)
    r <- matrix(c(z,z,z,bg), nrow=2)
  if(z >= 2 && z < 3)  # PC    (positive correlation)
    r <- matrix(c(bg,z,z,bg), nrow=2)
  if(z >= 3 && z < 4)  # LIL   (low x implies low y)
    r <- matrix(c(bg,z,z,z), nrow=2)
  if(z >= 4 && z < 5)  # UNL   (unspecified non-linear)
    r <- matrix(c(z,z,z,z), nrow=2)
  if(z >= 5 && z < 6)  # HIL   (high x implies low y)
    r <- matrix(c(z,z,bg,z), nrow=2)
  if(z >= 6 && z < 7)  # NC    (negative correlation)   
    r <- matrix(c(z,bg,bg,z), nrow=2)
  if(z >= 7 && z < 8)  # LIH   (low x implies low y)   
    r <- matrix(c(z,bg,z,z), nrow=2)
  r
}


# Extend this function to expand another matrix on the basis of another
expand.cls.2 <- function(CLS, TARGET=NULL, pad=FALSE, bg=NA) {
  if(is.null(TARGET))
    TARGET <- CLS
  if (pad)
    G <- mat.or.vec(nrow(CLS)*3+1, ncol(CLS)*3+1)
  else
    G <- mat.or.vec(nrow(CLS)*2, ncol(CLS)*2)
  G[,] <- bg # fill with background enumeration
  
  for(i in 0:(nrow(CLS)-1)) {
    for(j in 0:(ncol(CLS)-1)) {
      gly <- to.glyph(CLS[i+1,j+1], bg)
      if(is.na(bg))
        gly[!is.na(gly)] <- TARGET[i+1,j+1]
      else
        gly[gly!=bg] <- TARGET[i+1,j+1]
      if (pad)
        G[(i*3+2):(i*3+3),(j*3+2):(j*3+3)] <- gly
      else
        G[(i*2+1):(i*2+2),(j*2+1):(j*2+2)] <- gly
    }
  }
  G
}

GG <- expand.cls.2(CLS.ord, DELTA.ord, bg=NA)
GG.pad <- expand.cls.2(CLS.ord, DELTA.ord, bg=NA, pad=TRUE)
GG[1:4,1:4]
GG.pad[1:4,1:4]
heatmap.3(GG, MIN=min(DELTA), MAX=max(DELTA), reorder=FALSE)
## Show me preliminary recolored map
heatmap.3(GG.pad, MIN=min(DELTA), MAX=max(DELTA), reorder=FALSE)
R <- splom(CLS.ord, DCOR.ord)

## Look at corresponding splom


## Expand class of constant matrix, color solid
ZERO <- mat.or.vec(nrow(CLS), ncol(CLS))
GG.solid.pad <- expand.cls.2(CLS.ord, ZERO, bg=NA, pad=TRUE)

G <- GG.solid.pad
G[1:10,1:10]
Img <- t(G)[,seq(nrow(G),1,-1)]
w<-ncol(G); h<-nrow(G)
image(1:w, 1:h, Img, col=c("black"), breaks=c(-1,1), axes=FALSE, xlab="", ylab="")
draw.glyphs(G)


#### ------------------------------
## Recolor map
