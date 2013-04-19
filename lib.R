library("RColorBrewer")
library("gplots")

GRID.COL <-  "#ffffff"
CLS.ENUM <- list("0"="NA", "1"="HIH", "2"="PC", "3"="LIL", "4"="UNL", "5"="HIL", "6"="NC", "7"="LIH", "8"="PAD")
GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")
GLYPH.COLS.MAX <- c("#ffffff", "#00b271", "#0089d9", "#3424b3", "#000000", "#a3033c", "#d82a36", "#eb5218", "#ffffff")
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

# Append "" after each name in double-sized list
expand.names <- function(i, name.list) { 
  if (i%%2==1) {
    name.list[(i+1)/2]
  } else {
    ""
  }
}

# Fix enumeration as outputted by Python program
renumerate.fix <- function(BC) {
  BC[BC==0] <- -1
  BC[BC==7] <- 0
  BC[BC>=4 & BC<=6] <- BC[BC>=4 & BC<=6]+1
  BC[BC==-1] <- 4
  # fix 1/3 switch
  BC[BC==1] <- -1
  BC[BC==3] <- 1
  BC[BC==-1] <- 3
  BC
}

splom <- function(CLS, DCOR=NULL, asGlyphs=FALSE, pad=FALSE, grid="auto", grid.col=GRID.COL, lwd=1, ...) {
  if(is.null(DCOR)) {
    R <- splom.cls(CLS, asGlyphs=asGlyphs, pad=pad, ...)
  } else {
    R <- splom.dcor(CLS, DCOR=DCOR, asGlyphs=asGlyphs, pad=pad, ...)
  }
  if (grid==TRUE || (grid=="auto" && asGlyphs)) {
    if (asGlyphs)
      if (pad) 
        grid.offset <- 3
      else
        grid.offset <- 2
    else
      grid.offset <- 1
    w <- ncol(R$G); h <- nrow(R$G)
    draw.grid(grid.offset, w, h, grid.col, lwd)
  }
  R
}

# ========================================
# ========================================

make.color.bins <- function(N=15, high.sat=TRUE) {
  if(high.sat) {
    COLORS <- GLYPH.COLS.MAX
  } else {
    COLORS <- GLYPH.COLS
  }
  COLOR.M <- sapply(COLORS, function(color) colorRampPalette(c("#ffffff", color))(N+1))
  c(COLOR.M)
}
make.breaks <- function(MAX=1, MIN=0, N=15, MOST=1, LEAST=0) {
  th <- (MAX-MIN)/N   ## bin for above and below bin threshold
  c(LEAST, sapply(0:(N-1), function(i) MIN+i*th), MOST)
}
make.offsets <- function(breaks, N=15) {
  offsets <- sapply(1:(N+1), function(i) (breaks[i]+breaks[i+1])/2)
  offsets <- c(offsets, tail(offsets,1)) ## offset beyond max value is still max value
  offsets
}

get.order.cls.dcor <- function(CLS, DCOR, DCOR.weight=2, clust.meth="complete", CLS.enum.dist=F, DCOR.include.PCC=F) {
  R = list()
  if (CLS.enum.dist) {
    D.cls.r <- dist(CLS)
    D.cls.c <- dist(t(CLS))
  } else {
    D.cls.r <- gen.glyph.dist.m(CLS)
    D.cls.c <- gen.glyph.dist.m(t(CLS))
    # max dist: 4*m
  }
  if (DCOR.include.PCC) {
    D.DCOR.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)
    D.DCOR.c <- as.dist(1-cor(DCOR, method="pearson")) + dist(t(DCOR))
  } else {
    D.DCOR.r <- dist(DCOR)
    # max dist: sqrt(m)
    D.DCOR.c <- dist(t(DCOR))
  }

  Rowv <- rowMeans(DCOR, na.rm = TRUE)
  Colv <- colMeans(DCOR, na.rm = TRUE)
  R$D.row <- D.DCOR.r*DCOR.weight+sqrt(D.cls.r)
  R$D.col <- D.DCOR.c*DCOR.weight+sqrt(D.cls.c)
  R$Rhclust <- as.dendrogram(hclust(R$D.row, method=clust.meth))
  R$Rhclust <- reorder(R$Rhclust, Rowv)
  R$Chclust <- as.dendrogram(hclust(R$D.col, method=clust.meth))
  R$Chclust <- reorder(R$Chclust, Colv)
  R
}

## Draw only class enumerations.
splom.cls <- function(CLS, reorder=TRUE, asGlyphs=FALSE, pad=FALSE, ...) {
  R <- get.cls.order(CLS)
  rowInd <- order.dendrogram(R$Rhclust)
  colInd <- order.dendrogram(R$Chclust)
  
  w <- ncol(R$G); h <- nrow(R$G)
  sapply(1:(length(rownames(DCOR))*2), function(i) expand.names(i,rownames(DCOR)))
  
  if (reorder)
    R$G <- CLS[rowInd, colInd]
  else
    R$G <- CLS
  if (asGlyphs) 
    R$G <- expand.cls(R$G, pad=pad)
  draw.glyphs(R$G, grid=0, ...)
  R
}

## Draw DCOR scaled class enumerations.
splom.dcor <- function(CLS, DCOR, reorder=TRUE, asGlyphs=FALSE, pad=FALSE, N=15, MIN=0.1, MAX=0.8, MOST=1, LEAST=0, DCOR.weight=2, useRaster=FALSE, high.sat=TRUE, draw.labs=T, clust.meth="complete", ...) {
  ## Clustering
  if (reorder) {
    R <- get.order.cls.dcor(CLS, DCOR, DCOR.weight, clust.meth=clust.meth)
    rowInd <- order.dendrogram(R$Rhclust)
    colInd <- order.dendrogram(R$Chclust)
    CLS <- CLS[rowInd, colInd]
    DCOR <- DCOR[rowInd, colInd]
  } else {
    rowInd <- 1:dim(DCOR)[1]
    colInd <- 1:dim(DCOR)[2]
    R <- list()
  }
  ## Generate color/class bins
  R$COLOR.V <- make.color.bins(N, high.sat)
  R$breaks <- make.breaks(MAX, MIN, N, MOST, LEAST)
  R$offsets <- make.offsets(R$breaks, N)
  ## Select the greatest offset less than x.
  choose.offset <- function(x) R$offsets[tail(which(R$breaks<=x),1)]
  N.BREAKS <- c(sapply(0:8, function(x) rep(x,N+1)+R$breaks[1:(N+1)]), 9)

  ## Make image heatmap from matrices
  OFFSET <- apply(DCOR, c(1,2), choose.offset)
  R$G <- CLS+OFFSET
  if (asGlyphs) {
    R$G <- expand.cls(R$G, pad=pad)
    if (pad)
      R$G[R$G==8] <- 8.01
  }
  Img <- t(R$G)[,seq(nrow(R$G),1,-1)]

  # Row and column labels
  w <- ncol(R$G); h <- nrow(R$G)
  sapply(1:(length(rownames(DCOR))*2), function(i) expand.names(i,rownames(DCOR)))

  if (asGlyphs) {
    labRow <- sapply(1:(length(rownames(DCOR))*2), function(i) expand.names(i,rev(rownames(DCOR))))
    labCol <- sapply(1:(length(colnames(DCOR))*2), function(i) expand.names(i,colnames(DCOR)))
  } else {
    labRow <- rev(rownames(DCOR))
    labCol <- colnames(DCOR)
  }

  ## Draw image
  nc <- dim(DCOR)[2]
  nr <- dim(DCOR)[1]
  if (draw.labs)
    par(mar = c(5, 0, 5, 5))
  image(1:w, 1:h, Img, xlab="", ylab="", col=R$COLOR.V, breaks=N.BREAKS, axes=FALSE, useRaster=useRaster, ...)
  if (draw.labs) {
    if (asGlyphs) {
      axis(1, 1:(nc*2), labels = labCol, las = 2, line = -0.5, tick = 0, 
           cex.axis = 0.7)
      axis(4, 1:(nr*2), labels = labRow, las = 2, line = -0.5, tick = 0, 
           cex.axis = 0.7)
    } else {
      axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
           cex.axis = 0.7)
      axis(4, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, 
           cex.axis = 0.7)
    }
  }
  
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(DCOR))) 
      (1:nr)[rowInd]
    else rownames(DCOR)
  
  R
}


## Extend heatmap.2 to generate a heatmap with decent colors and clustering.
## --------------------
## EXAMPLE:
## H <- heatmap.3(BC0.Dcor, ylab="CpG", xlab="mRNA")
heatmap.3 <- function(M, MIN=0.08, MAX=0.8, cols=brewer.pal(8,"RdYlBu"), reorder=TRUE, symm=FALSE, symkey=FALSE, key=TRUE, trace="none", clust.meth="average", ...) {
  heatmap_breaks <- seq(MIN,MAX,0.01)
  heatmap_cols <- rev(colorRampPalette(cols)(length(heatmap_breaks)-1))
  if (reorder) {
    # distance as sum of euclidean and correlation. WARNING: cor on columns, dist on rows!
    D.r <- as.dist(1-cor(t(M), method="pearson")) + dist(M)
    D.c <- as.dist(1-cor(M, method="pearson")) + dist(t(M))
    # Row and column ordering based on mean values
    Rowv <- rowMeans(M, na.rm = TRUE)
    Colv <- colMeans(M, na.rm = TRUE)
    Rhclust <- as.dendrogram(hclust(D.r, method=clust.meth))
    Rhclust <- reorder(Rhclust, Rowv)
    Chclust <- as.dendrogram(hclust(D.c, method=clust.meth))
    Chclust <- reorder(Chclust, Colv)
    dendrogram <- "both"
  } else {
    Rhclust <- NULL; Chclust <- NULL;
    dendrogram <- "none"
  }

  heatmap.2(as.matrix(M),
    col=heatmap_cols, breaks=heatmap_breaks,
    key=key, symkey=symkey, trace=trace, dendrogram=dendrogram,
    Rowv=Rhclust, Colv=Chclust, ...
  )
}


## Cluster an enumerated boolean dependency class matrix.
## --------------------
# Helper function to return class ordering as Row / Column dendrogram objs.
get.cls.order <- function(CLS, clust.meth="complete") {
  R = list()
  D.cls.r <- gen.glyph.dist.m(CLS)
  D.cls.c <- gen.glyph.dist.m(t(CLS))
  Rowv <- rowMeans(CLS, na.rm = TRUE)
  Colv <- colMeans(CLS, na.rm = TRUE)
  R$Rhclust <- as.dendrogram(hclust(D.cls.r, method=clust.meth))
  R$Rhclust <- reorder(R$Rhclust, Rowv)
  R$Chclust <- as.dendrogram(hclust(D.cls.c, method=clust.meth))
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


## Given a class matrix, construct a glyph matrix.
## --------------------
expand.cls <- function(CLS, TARGET=NULL, pad=FALSE, bg=NA) {
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

# 
summary.plots <- function(CLS, DCOR, sym=F) {
  if (sym) {
    CLS <- CLS[upper.tri(CLS, diag = FALSE)]
    DCOR <- DCOR[upper.tri(DCOR, diag = FALSE)]
    upper.note <- " (upper triangle)"
  } else {
    upper.note <- ""
  }
  ENUM <- list()
  for (i in 0:7) ENUM[[as.character(i)]] <- 0
  Z <- split(DCOR, CLS)
  for (n in names(Z)) ENUM[[n]] <- Z[[n]]
  names(ENUM) <- CLS.ENUM[match(names(ENUM), names(CLS.ENUM))]
  boxplot(ENUM, col=GLYPH.COLS[1:8], main=paste0("dCOR per Class", upper.note))
  barplot(sapply(ENUM,length), col=GLYPH.COLS[1:8], main=paste0("Boolean Class Frequency", upper.note))
  hist(as.matrix(DCOR), main=paste0("Histogram of all-pairs dCOR", upper.note))
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
  draw.grid(grid, w, h, grid.col, lwd)
}

draw.grid <- function(grid, w, h, grid.col=GRID.COL, lwd=1) {
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

## TODO: Manually Generate intensity scale
## -----

# Load distances
#G.DIST.TABLE <- read.table("glyph.dists.csv", as.is=T, sep=",")
cls.enum.names <- c("hih","pc","lil","unl","hil","nc","lih","na")
r1 <- c(0,1,2,1,2,3,2,2)
r2 <- c(1,0,1,2,3,4,3,2)
r3 <- c(2,1,0,1,2,3,2,2)
r4 <- c(1,2,1,0,1,2,1,2)
r5 <- c(2,3,2,1,0,1,2,2)
r6 <- c(3,4,3,2,1,0,1,2)
r7 <- c(2,3,2,1,2,1,0,2)
r8 <- c(2,2,2,2,2,2,2,0)
G.DIST.TABLE <- rbind(r1,r2,r3,r4,r5,r6,r7,r8)
rownames(G.DIST.TABLE) <- cls.enum.names
colnames(G.DIST.TABLE) <- cls.enum.names


glyph.dist.f <- function(A,B) {
  # Requires that NA is mapped to 8 rather than 0
  sum(apply(cbind(A,B), 1, function(p) G.DIST.TABLE[p[1],p[2]]))
}

# all-rows sum glyph hamming distance matrix
gen.glyph.dist.m <- function(BC, recast.na.0=T) {
  if(recast.na.0)
    BC[BC==0] <- 8 # R indexes from 1, not 0, so remap 0 to 1+ the biggest glyph enum (8=7+1)
  n <- nrow(BC)
  as.dist(outer(1:n,1:n, FUN = Vectorize( function(i,j) glyph.dist.f(BC[i,],BC[j,]) )))
}

