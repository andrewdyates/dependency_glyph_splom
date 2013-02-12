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

#GLYPH.COLS.MAX.2 <- c("#ffffff", "#00b371", "#0089d9", "#3424b3", "#000000", "#bf0063", "#bf000d", "#f04400", "#ffffff")
GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")
GLYPH.COLS.MAX <- c("#ffffff", "#00b271", "#0089d9", "#3424b3", "#000000", "#a3033c", "#d82a36", "#eb5218", "#ffffff")
#9a0049 ## more blue

# Row is level (1 is white, N is color); Col is glyph in GLYPH.COLS order


# There are N equally-spaced tints from MIN to MAX.
#   tints do not include white but do including the untinted color.
#   values under MIN are assigned the background color (white)
#   The 1st tint histogram bin (over MIN) is a tint; it is not white.
#   The last histogram bin is to MAX. Values over MAX are in this bin.
# values in the first bin are 
N <- 15 # number of tints including the color itself. Must be >= 1.
## First interpolation level is white for values below MIN.
## Thus there are N-1 tints, plus white, for each range.
## We assign values of at least MAX the untinted color.
## --------------------
## COLOR.V assigns the 1st N elements to GLYPH.COLS[1] tints,
##   the 2nd N elements to GLYPH.COLS[2] tints... etc.
##   thus, each color "column" gets N+1 entries: white, 9 tints, and the color
#COLOR.M <- sapply(GLYPH.COLS, function(color) colorRampPalette(c("#ffffff", color))(N+1))
COLOR.M <- sapply(GLYPH.COLS.MAX, function(color) colorRampPalette(c("#ffffff", color))(N+1))
# may adjust this to improve saturation
COLOR.V <- c(COLOR.M)

MIN <- 0.08
MAX <- max(DCOR) # 0.792174

# Mark class ordinals using a decimal (less than 1). Order by increasing dep strength.

## "Lower edge" of color bins. Bin is [i,i+1) of `breaks`, indexed from 1.
##    "0" goes into the first bin: [-0.5, 0.5)
##    "1" goes into second bin:    [0.5, 1.5)
##    ...
##    "7" goes into ninth bin (eight bins, plus lowest bound)  [7.5, 8.5)

th <- (MAX-MIN)/N   ## bin for above and below bin threshold
# do not add a break at MAX; extend that bin to end at the global max (1)
breaks <- c(0, sapply(0:(N-1), function(i) MIN+i*th), 1)

## Select greatest break less than x

offsets <- sapply(1:(N+1), function(i) (breaks[i]+breaks[i+1])/2)
offsets <- c(offsets, tail(offsets,1)) ## offset beyond max value is still max value
## Select the greatest offset less than x.
f <- function(x) offsets[tail(which(breaks<=x),1)]

## Generate offset fraction from DCOR histogram.
OFFSET <- apply(DCOR.ord, c(1,2), f)


G <- CLS.ord+OFFSET
Img <- t(G)[,seq(nrow(G),1,-1)]
# 
# convert COLOR.M to vector
# draw.glyphs(CLS.ord)
## ------------------------------
# generate color level breaks from 0:7 with decimal histograms

N.BREAKS <- c(sapply(0:8, function(x) rep(x,N+1)+breaks[1:(N+1)]), 9)

w<-ncol(G); h<-nrow(G)
image(1:w, 1:h, Img, col=COLOR.V, breaks=N.BREAKS, axes=FALSE, xlab="", ylab="", useRaster=FALSE)



## ==============================
## Cluster on offset class enumerations
## ==============================
R <- get.cls.order(G)
rowInd <- order.dendrogram(R$Rhclust)
colInd <- order.dendrogram(R$Chclust)
G.ord <- G[rowInd, colInd]
G.CLS.ord <- CLS.ord[rowInd, colInd]
G.DCOR.ord <- DCOR.ord[rowInd, colInd]

Img2 <- t(G.ord)[,seq(nrow(G.ord),1,-1)]
image(1:w, 1:h, Img2, col=COLOR.V, breaks=N.BREAKS, axes=FALSE, xlab="", ylab="", useRaster=FALSE)

## ==============================
## Cluster on class enumerations, sort by dCOR
## ==============================
CLS  <- BC0.cls
DCOR <- BC0.dcor

D.cls.r <- dist(CLS)
D.cls.c <- dist(t(CLS))
Rowv <- rowMeans(DCOR, na.rm = TRUE)
Colv <- colMeans(DCOR, na.rm = TRUE)
R$Rhclust <- as.dendrogram(hclust(D.cls.r, method="average"))
R$Rhclust <- reorder(R$Rhclust, Rowv)
R$Chclust <- as.dendrogram(hclust(D.cls.c, method="average"))
R$Chclust <- reorder(R$Chclust, Colv)

rowInd <- order.dendrogram(R$Rhclust)
colInd <- order.dendrogram(R$Chclust)
CLS.ord <- CLS[rowInd, colInd]
DCOR.ord <- DCOR[rowInd, colInd]

OFFSET <- apply(DCOR.ord, c(1,2), f)
G <- CLS.ord+OFFSET
Img <- t(G)[,seq(nrow(G),1,-1)]
image(1:w, 1:h, Img, col=COLOR.V, breaks=N.BREAKS, axes=FALSE, xlab="", ylab="", useRaster=FALSE)


## ==============================
## Cluster on class enumerations w/ dCOR dist
##   sort by dCOR
## ==============================
CLS  <- BC0.cls
DCOR <- BC0.dcor
DCOR.W <- 2 # best by rough experiment

D.cls.r <- dist(CLS)
D.cls.c <- dist(t(CLS))
D.DCOR.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)
D.DCOR.c <- as.dist(1-cor(DCOR, method="pearson")) + dist(t(DCOR))

Rowv <- rowMeans(DCOR, na.rm = TRUE)
Colv <- colMeans(DCOR, na.rm = TRUE)
Rhclust <- as.dendrogram(hclust(D.DCOR.r*DCOR.W+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
Chclust <- as.dendrogram(hclust(D.DCOR.c*DCOR.W+D.cls.c, method="average"))
Chclust <- reorder(Chclust, Colv)

rowInd <- order.dendrogram(Rhclust)
colInd <- order.dendrogram(Chclust)
CLS.ord <- CLS[rowInd, colInd]
DCOR.ord <- DCOR[rowInd, colInd]

OFFSET <- apply(DCOR.ord, c(1,2), f)
G <- CLS.ord+OFFSET
Img <- t(G)[,seq(nrow(G),1,-1)]
image(1:w, 1:h, Img, col=COLOR.V, breaks=N.BREAKS, axes=FALSE, xlab="", ylab="", useRaster=FALSE)


### Get list of dependencies by dependency category
### ------------------------------
# Load mRNA
methyl <- as.matrix(read.table("data/Methyl_correct_aligned.tab.gz", header=TRUE, sep="\t", row.names=1))
mRNA <- as.matrix(read.table("data/mRNA_correct_aligned.tab.gz", header=TRUE, sep="\t", row.names=1))
# load("data/methyl_mrna.RData")

#CLS  <- BC0.cls
#DCOR <- BC0.dcor
CLS  <- BCBig.cls
DCOR <- BCBig.dcor

# Get methyl and mRNA rows included in CLS and DCOR aligned to rows and columns
#   rows: CpG
#   cols: mRNA
ROW.VARS = methyl[which(rownames(methyl) %in% rownames(CLS)),]
COL.VARS = mRNA[which(rownames(mRNA) %in% colnames(CLS)),]

# Plot distribution of class and dCOR
# ------------------------------
boxplot(split(DCOR, CLS), col=GLYPH.COLS[2:8])
barplot(sapply(split(DCOR, CLS),length), col=GLYPH.COLS[2:8])
hist(DCOR)




DCOR[CLS == 4]
qq <- which(CLS == 4, arr.ind=T)
ww <- order(DCOR[qq], decreasing=TRUE)
## qq[ww[1],]
## row col 
##  76  16
DCOR[qq[ww[1],]]

DCOR[rbind(qq[ww[1],])]
## !!NOTE: DCOR values in Python are energy dCOR values squared.

# color dots by tissue
col <- rep("black", length(X)) # FCTX
col[grep("TCTX", names(X))] <- "red"
col[grep("CRBLM", names(X))] <- "blue"
col[grep("PONS", names(X))] <- "green"

test <- function(n) {
for(i in 1:n) {
  x <- qq[ww[i],][1]
  y <- qq[ww[i],][2]
  X <- ROW.VARS[x,]
  Y <- COL.VARS[y,]
  png(paste0("nonlin",i,".",x,".",y,".png"))
  plot(X,Y, col=col)
  dev.off()
  print(paste0(cor(X,Y), ": ", dcor(X,Y)))
}
}



test2 <- function(n) {
for(cls in c(1,2,3,4,5,6,7)) {
  cls.idx <- which(CLS == cls, arr.ind=T)
  dcor.rank <- order(DCOR[cls.idx], decreasing=TRUE)
  print(cls)
  for(i in 1:n) {
    x <- cls.idx[dcor.rank[i],][1]
    y <- cls.idx[dcor.rank[i],][2]
    X <- ROW.VARS[x,]
    Y <- COL.VARS[y,]
    pcc <- cor(X,Y)
    dcc <- dcor(X,Y)
    png(paste0("cls.",cls,".",i,".",x,".",y,"pcc",pcc,"dcor",dcc,".png"))
    plot(X,Y, col=col)
    dev.off()
    print(paste0(pcc, ": ", dcc))
  }
}
}
