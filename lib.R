library("RColorBrewer")
library("gplots")
#library("flashClust")
library("fastcluster")

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

get.order.cls.dcor <- function(CLS, DCOR, DCOR.weight=2, clust.meth="average", CLS.enum.dist=F, DCOR.include.PCC=F, sym=F, col.cls.dist=NULL, row.cls.dist=NULL) {
  R = list()
  if (CLS.enum.dist) {
    D.cls.r <- dist(CLS)
    if (sym)
      D.cls.c <- D.cls.r
    else
      D.cls.c <- dist(t(CLS))
  } else {
    D.cls.r <- gen.glyph.dist.m(CLS, precomp=row.cls.dist)
    if (sym)
      D.cls.c <- D.cls.r
    else
      D.cls.c <- gen.glyph.dist.m(t(CLS), precomp=col.cls.dist)
    # max dist: 4*m
  }
  if (DCOR.include.PCC) {
    D.DCOR.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)
    if (sym)
      D.DCOR.c <- D.DCOR.r
    else
      D.DCOR.c <- as.dist(1-cor(DCOR, method="pearson")) + dist(t(DCOR))
  } else {
    D.DCOR.r <- dist(DCOR)
    if (sym)
      D.DCOR.c <- D.DCOR.r
    else
      D.DCOR.c <- dist(t(DCOR))
    # max dist: sqrt(m)
  }

  Rowv <- rowMeans(DCOR, na.rm = TRUE)
  R$D.row <- D.DCOR.r*DCOR.weight+sqrt(D.cls.r)
  R$Rhclust <- as.dendrogram(hclust(R$D.row, method=clust.meth))
  R$Rhclust <- reorder(R$Rhclust, Rowv)

  if (sym) {
    Colv <- Rowv
    R$D.col <- R$D.row
    R$Chclust <- R$Rhclust
  } else {
    Colv <- colMeans(DCOR, na.rm = TRUE)
    R$D.col <- D.DCOR.c*DCOR.weight+sqrt(D.cls.c)
    R$Chclust <- as.dendrogram(hclust(R$D.col, method=clust.meth))
    R$Chclust <- reorder(R$Chclust, Colv)
  }
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
splom.dcor <- function(CLS, DCOR, reorder=TRUE, asGlyphs=FALSE, pad=FALSE, N=15, MIN=0.2, MAX=1, MOST=1, LEAST=0, DCOR.weight=2, useRaster=FALSE, high.sat=TRUE, draw.labs=T, clust.meth="average", sym=F, row.cls.dist=NULL, col.cls.dist=NULL, ...) {
  ## Clustering
  if (reorder) {
    R <- get.order.cls.dcor(CLS, DCOR, DCOR.weight, clust.meth=clust.meth, sym=sym, row.cls.dist=row.cls.dist, col.cls.dist=col.cls.dist)
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
get.cls.order <- function(CLS, clust.meth="average") {
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
  summary.plots.vector(CLS, DCOR)
}

get.size <- function(s) {
  if (length(s)==1 && is.na(s)) {
    0
  } else {
    length(s)
  }
}
summary.plots.vector <- function(CLS, DCOR) {
  ENUM <- list()
  for (i in 0:7) ENUM[[as.character(i)]] <- NA
  Z <- split(DCOR, CLS)
  for (n in names(Z)) ENUM[[n]] <- Z[[n]]
  names(ENUM) <- CLS.ENUM[match(names(ENUM), names(CLS.ENUM))]
  if (length(DCOR) < 1000000) {
    boxplot(ENUM, col=GLYPH.COLS[1:8], main="dCOR per Class", ylim=c(0,1))
  } else {
    ENUM.samp <- list()
    qq <- sample.int(length(DCOR), size=1000000)
    Z.samp <- split(DCOR[qq], CLS[qq])
    for (n in names(Z.samp)) ENUM.samp[[n]] <- Z.samp[[n]]
    names(ENUM.samp) <- CLS.ENUM[match(names(ENUM.samp), names(CLS.ENUM))]
    boxplot(ENUM.samp, col=GLYPH.COLS[1:8], main="dCOR per Class (1e6 sample)", ylim=c(0,1))
  }
  sizes <- sapply(ENUM,get.size)
  barplot(sizes, col=GLYPH.COLS[1:8], main="Boolean Class Frequency")
  hist(as.matrix(DCOR), main="Histogram of all-pairs dCOR", xlim=c(0,1), breaks=0:10/10)
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
make.cls.dist.M <- function() {
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
  G.DIST.TABLE
}
G.DIST.TABLE <- make.cls.dist.M()


glyph.dist.f <- function(A,B) {
  # Requires that NA is mapped to 8 rather than 0
  sum(apply(cbind(A,B), 1, function(p) G.DIST.TABLE[p[1],p[2]]))
}

# all-rows sum glyph hamming distance matrix
# WARNING: THIS IS EXTREMELY SLOW
gen.glyph.dist.m <- function(BC, recast.na.0=T, precomp=NULL) {
  if(!is.null(precomp)) {
    stopifnot(dim(precomp)[1]==dim(BC)[1])
    stopifnot(dim(precomp)[2]==dim(BC)[1])
    return(precomp);
  } else {
    if (class(BC) != "matrix") BC <- as.matrix(BC)
    if(recast.na.0)
      BC[BC==0] <- 8 # R indexes from 1, not 0, so remap 0 to 1+ the biggest glyph enum (8=7+1)
    n <- nrow(BC)
    as.dist(outer(1:n,1:n, FUN = Vectorize( function(i,j) glyph.dist.f(BC[i,],BC[j,]) )))
  }
}




### --------------------------------------------------
### COHERENCE COLLAPSE

# coherences is max 1-(1/4) mean hamming dist to single class
# BOOL_ENUM = {0:'NA', 1:'XiY', 2:'PC', 3:'YiX', 4:'UNL', 5:'MX', 6:'NC', 7:'OR'}
# NOTE: NA:8 in R CLS dist matrix
# Count occurences of classes. Indexed from 1. NA:0 is index 8.
get.cls.counts <- function(CLS, sym=T) {
  if (sym) {
    stopifnot(dim(CLS)[1] == dim(CLS)[2])
    C <- CLS[upper.tri(CLS)]
  } else {
    C <- c(CLS)
  }
  C <- na.omit(c(C))
  # indexed from 1. NA is index 8
  if(any(C==0)) stop("Set NA CLS==0 to 8.")
  cls.counts <- sapply(c(1:8), function(i) sum(C==i))
  # account for arbitrary symmetry across diagonal for XiY and YiX classes
  if (sym) {
    cls.counts[1] = cls.counts[1] + cls.counts[3]
    cls.counts[3] = 0
  }
  cls.counts
}
# Calculate coherence for a cls given a count of other classes.
#   cls.counts enumerated from 1. NA is index 8 (not index 0)
get.coh <- function(cls.counts, cls) {
  stopifnot(length(cls.counts) == 8)
  stopifnot(cls >= 1 && cls <= 8)
  n <- sum(cls.counts)
  sigma <- sum(cls.counts * G.DIST.TABLE[cls,])
  e <- 1-(sigma/n/4)
  (e-0.5)*2 # center at 0, range -1 to 1
}
# Calculate coherence coherence vector for all classes
#   cls.counts enumerated from 1. NA is index 8 (not index 0)
get.coh.vec <- function(CLS, sym=T) {
  if (length(CLS) == 1) sym<-F
  cls.counts <- get.cls.counts(CLS, sym)
  sapply(1:8, function(i) get.coh(cls.counts, i))
}

# Return most coherent class enumeration given coherence vector.
#   coh.v enumerated from 1. NA is index 8 (not index 0)
choose.coh.cls <- function(coh.v, min.coh=0.74) {
  max.coh <- max(coh.v)
  max.coh.i <- coh.v == max.coh
  R <- list(); R$coh <- max.coh
  # maximum coherence below threshold. Return UNL
  b <- max.coh < min.coh
  if (is.na(b)) {
    print(coh.v)
    stop("BAD")
  }
  if (b) {
    R$cls <- 4
  } else {
    if (sum(max.coh.i)==1) {
      R$cls <- which(max.coh.i)
    # Handle ties.
    # 1. ignore NA class
    } else if (sum(max.coh.i[1:7])==1) {
      R$cls <- which(max.coh.i[1:7])
    # 1.5 ignore NA and UNL class
    } else if (sum(max.coh.i[c(1,2,3,5,6,7)]) == 1) {
      max.coh.i[4] <- FALSE
      R$cls <- which(max.coh.i[1:7])
    # 2. if ties for opposite signs, return UNL
    } else if (any(max.coh.i[c(1,2,3)]) & any(max.coh.i[c(5,6,7)])) {
      R$cls <- 4
    # 3. if all of a sign or both asyms of a sign, return sym of that sign
    } else if (all(max.coh.i[c(1,2,3)]) | all(max.coh.i[c(1,3)])) {
      R$cls <- 2
    } else if (all(max.coh.i[c(5,6,7)]) | all(max.coh.i[c(5,7)])) {
      R$cls <- 6
    # 4. if an asym and sym of the same sign, return asym
    } else if (all(max.coh.i[c(1,2)]) | all(max.coh.i[c(2,3)])) {
      R$cls <- 2
    } else if (all(max.coh.i[c(5,6)]) | all(max.coh.i[c(6,7)])) {
      R$cls <- 6
    }
  }
  R
}

get.mean.dcor <- function(DCOR, sym=F) {
  if (sym) {
    D <- DCOR[upper.tri(DCOR)]
  } else {
    D <- DCOR
  }
  u<-mean(D)
  if(is.na(u))
    u <- DCOR # mean is value itself
  u
}

# Return collapsed class, coherence, and overall coherence by indices.
collapse.cls <- function(CLS, idx, DCOR=NULL, dcor.sig=NULL) {
  CLS[CLS==0] <- 8
  idx <- as.factor(idx)
  n <- length(levels(idx))
  R <- list()
  R$idx <- idx
  R$SIZE <- matrix(0, nrow=n, ncol=n)
  R$CLS <- matrix(0, nrow=n, ncol=n)
  rownames(R$CLS) <- 1:n
  colnames(R$CLS) <- 1:n
  R$COH <- matrix(0, nrow=n, ncol=n)
  R$MIX.SIGN <- matrix(0.0, nrow=n, ncol=n)
  R$MIX.DIR <- matrix(0.0, nrow=n, ncol=n)
  R$LOSER.CLUST.EDGES <- rep(0,n)
  R$members <- split(rownames(CLS), idx)
  if (!is.null(DCOR)) {
    R$DCOR <- matrix(0, nrow=n, ncol=n)
    rownames(R$DCOR) <- 1:n
    colnames(R$DCOR) <- 1:n
    R$DCOR.mid <- (dcor.sig*2 + 1)/3
  } else {
    R$DCOR <- NULL
  }
  # compute overall coherence
  for (i in 1:length(levels(idx))) {
    for (j in 1:length(levels(idx))) {
      iv <- idx == levels(idx)[i]
      jv <- idx == levels(idx)[j]
      C <- CLS[iv,jv]
      if (!is.null(DCOR) && !is.null(dcor.sig)) {
        C[DCOR[iv,jv]<dcor.sig] <- NA
      }
      if (!is.null(DCOR)) {
        R$DCOR[i,j] <- get.mean.dcor(DCOR[iv,jv], sym=i==j)
        # if this is a cluster, does it contain a weakly-dependent or negative class "loser" feature?
        if (i==j)
          R$LOSER.CLUST.EDGES[i] <- sum((DCOR[iv,jv] < R$DCOR.mid & C!=2) | (C %in% c(5,6,7)) | DCOR[iv,jv] < dcor.sig) / 2
      } else {
        if (i==j)
          R$LOSER.CLUST.EDGES[i] <- sum(C %in% c(5,6,7))
      }
      if (!is.null(DCOR) && !is.null(dcor.sig) && R$DCOR[i,j] < dcor.sig) {
        R$MIX.SIGN[i,j] <- 0
        R$MIX.DIR[i,j]  <- 0
        R$CLS[i,j] <- 4
        R$COH[i,j] <- 1
      } else {
        r <- choose.coh.cls(get.coh.vec(C, sym=i==j))
        if (is.null(r$cls)) {
          print(C)
          print(r)
          stop("Null class.")
        }
        R$CLS[i,j] <- r$cls
        R$COH[i,j] <- r$coh
        R$MIX.SIGN[i,j] <- count.mix.sign(na.omit(c(C)))
        R$MIX.DIR[i,j]  <- count.mix.pos.dir(na.omit(c(C))) + count.mix.neg.dir(na.omit(c(C)))
      }
      if (i == j) {
        stopifnot(sum(iv)==sum(jv))
        R$SIZE[i,j] <- sum(iv)*(sum(iv)-1)/2
      } else {
        R$SIZE[i,j] <- sum(iv)*sum(jv)
      }
    }
  }
  # get matrix of "critical flaws"
  SUM.FLAWS <- R$MIX.SIGN + R$MIX.DIR
  R$CRIT <- (log2(R$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
  R
}

# count number of cross edges
count.mix.sign <- function(CLS)
  min(sum(CLS %in% c(1,2,3)), sum(CLS %in% c(5,6,7)))

count.mix.pos.dir <- function(CLS)
  min(sum(CLS==1), sum(CLS==3))

count.mix.neg.dir <- function(CLS)
  min(sum(CLS==5), sum(CLS==7))

has.mix.sign <- function(CLS)
  any(c(1,2,3) %in% CLS) & any(c(5,6,7) %in% CLS)

has.mix.direction <- function(CLS)
  all(c(1,3) %in% CLS) | all(c(5,7) %in% CLS)

get.wavg.score <- function(coh, size) {
  nn <- sum(size)
  p <- size/nn
  sum(p*coh)
}

get.coh.M.score <- function(COLLAPSED, min.dcor=0) {
  R <- list()
  all.tri <- upper.tri(COLLAPSED$COH,diag=F)
  all.diag <- upper.tri(COLLAPSED$COH,diag=T) & !all.tri
  if (is.null(COLLAPSED$DCOR)) {
    tri <- all.tri
    diag <- all.diag
    R$dcor.clust <- NULL
  } else {
    d <- COLLAPSED$DCOR >= min.dcor
    tri <- d & all.tri
    diag <- d & all.diag
    R$min.dcor <- min.dcor
    R$num.dcor.clust <- sum(d) # significant edges 
  }
  R$all.wavg <- get.wavg.score(COLLAPSED$COH, COLLAPSED$SIZE)
  R$tri.wavg <- get.wavg.score(COLLAPSED$COH[tri], COLLAPSED$SIZE[tri])
  R$all.tri.wavg <- get.wavg.score(COLLAPSED$COH[all.tri], COLLAPSED$SIZE[all.tri])
  R$diag.wavg <- get.wavg.score(COLLAPSED$COH[diag], COLLAPSED$SIZE[diag])
  R$all.diag.wavg <- get.wavg.score(COLLAPSED$COH[all.diag], COLLAPSED$SIZE[all.diag])
  R$sum.edge.xsign <- sum(COLLAPSED$MIX.SIGN[tri]) # sum total of sign flaws
  R$sum.edge.xdir <- sum(COLLAPSED$MIX.DIR[tri]) # sum total of xdir flaws
  R$sum.edge.flaws <- R$sum.edge.xsign + R$sum.edge.xdir
  z <- (COLLAPSED$MIX.SIGN[tri]>0) | (COLLAPSED$MIX.DIR[tri]>0)
  R$edge.flaws <- sum(z)
  zz <- (COLLAPSED$MIX.SIGN[tri]>0)
  R$edge.sign.flaws <- sum(zz)
  zz <- (COLLAPSED$MIX.DIR[tri]>0)
  R$edge.dir.flaws <- sum(zz)
  R$edge.all.n <- sum(all.tri)
  R$edge.n <- sum(tri)
  R$clust.all.n <- sum(all.diag)
  R$clust.n <- sum(diag)
  R$loser.clusters <- sum(COLLAPSED$LOSER.CLUST.EDGES>1)
  R
}

# ignore contribution of edges below dCor threshold
get.compression <- function(CLS, H, DCOR=NULL, min.dcor=0, max.k=NULL) {
  if (is.null(max.k)) 
    n <- round(dim(CLS)[1]/2)+2
  else
    n <- max.k
  R <- list()
  #R$clust.wavg <- rep(0,n)
  #R$edge.wavg <- rep(0,n)
  #R$edge.flaws <- rep(0,n)
  for(k in 1:n) {
    idx <- cutree(H,k)
    C <- collapse.cls(CLS,idx,DCOR)
    S <- get.coh.M.score(C, min.dcor)
    R$all.wavg[k] <- S$all.wavg
    R$clust.wavg[k] <- S$diag.wavg
    R$edge.wavg[k] <- S$tri.wavg
    
    R$all.clust.wavg[k] <- S$all.diag.wavg
    R$all.edge.wavg[k] <- S$all.tri.wavg

    R$edge.flaws[k] <- S$edge.flaws
    R$edge.sign.flaws[k] <- S$edge.sign.flaws
    R$edge.sum.flaws[k] <- S$edge.sum.flaws
    R$edge.xsign[k] <- S$edge.xsign
    R$edge.xdir[k] <- S$edge.xdir
    R$edge.n[k] <- S$edge.n
    R$edge.all.n[k] <- S$edge.all.n
    R$clust.n[k] <- S$clust.n
    R$clust.all.n[k] <- S$clust.all.n
    if (!is.null(S$dcor.clust))
      R$dcor.clust[k] <- S$dcor.clust
    # REPORT
    cat("k: ", k, "\n")
    cat("tri.wavg: ", S$tri.wavg, "\n")
    cat("extant edge sign-only flaws: ",S$edge.sign.flaws, "\n")
    cat("extant edge flaws: ",S$edge.flaws, "\n")    
    cat("sum edge flaws: ",S$edge.sum.flaws, "\n")
    cat("edge.n: ", S$edge.n, "\n")
    cat("edge.all.n: ", S$edge.all.n, "\n")
    cat("S$edge.flaws / S$edge.all.n: ", S$edge.flaws / S$edge.all.n, "\n")
    cat("S$edge.flaws / S$edge.n: ", S$edge.flaws / S$edge.n, "\n")
    print("---")
  }
  R
}

## --------------------------------------------------
## COLLAPSE WEAKS
## --------------------------------------------------
# 0: no class; 1: and; 2: rn4c (row necessary for col); 3: cn4r (col necessary for row); 4: xor; 5: mix, 6: no class

collapse.weak <- function(WEAK, idx) {
  idx <- as.factor(idx)
  n <- length(levels(idx))
  COMP.W <- matrix(0, nrow=n, ncol=n)
  rownames(COMP.W) <- 1:n
  colnames(COMP.W) <- 1:n
  for (i in 1:length(levels(idx))) {
    for (j in 1:length(levels(idx))) {
      iv <- idx == levels(idx)[i]
      jv <- idx == levels(idx)[j]
      W <- WEAK[iv,jv]
      w <- choose.weak(W, sym=i==j)
      COMP.W[i,j] <- w
    }
  }
  COMP.W
}

choose.weak <- function(WEAK, sym=F) {
  cnt <- weak.counts(WEAK, sym=sym)
  n <- sum(cnt)
  max.cnt <- max(cnt)
  max.cnt.i <- cnt == max.cnt
  if (max.cnt < n/2) {
    return(5)
  } else if (max.cnt.i[2]) {
    if (cnt[3] == 0 && cnt[4] == 0)
      return(2)
    else
      return(5)
  } else if (max.cnt.i[3]) {
    if (cnt[2] == 0 && cnt[4] == 0)
      return(3)
    else
      return(5)
  } else if (max.cnt.i[4]) {
    if(cnt[2] == 0 && cnt[3] == 0 && cnt[1] == 0)
      return(4)
    else
      return(5)
  } else if (max.cnt.i[1]) {
    if(cnt[4] == 0)
      return(1)
    else
      return(5)
  } else {
    return(5)
  }
  print("WARNING: edge case in choose.weak.")
  5
}
      
weak.counts <- function(WEAK, sym=F) {
  WEAK[WEAK==0] <- 6
  if(sym) {
    stopifnot(dim(WEAK)[1] == dim(WEAK)[2])
    W <- WEAK[upper.tri(WEAK)]
  } else {
    W <- c(WEAK)
  }
  weak.counts <- sapply(c(1:6), function(i) sum(W==i))
  # account for arbitrary symmetry across diagonal
  if (sym) {
    weak.counts[2] = weak.counts[2] + weak.counts[3]
    weak.counts[3] = 0
  }
  weak.counts
}


