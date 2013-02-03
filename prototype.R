## 
library("RColorBrewer")
library(modeest)
library("gplots")

BC0.cls <- read.table("data/bc0.cls.tab", header=TRUE, sep=" ", row.names=1);
BC0.dcor <- read.table("data/bc0.dcor.tab", header=TRUE, sep=" ", row.names=1);
BCBig.cls <- read.table("data/bcBig.cls.tab.gz", header=TRUE, sep=" ", row.names=1);
BCBig.dcor <- read.table("data/bcbig.dcor.tab.gz", header=TRUE, sep=" ", row.names=1);
# ---------------
# Re-enumerate glyph symbols to put ND in the center, 0 as bg
renumerate <- function(BC) {
  BC[BC==0] <- -1
  BC[BC>=1 & BC<=3] <- BC[BC>=1 & BC<=3] -1
  BC[BC==-1] <- 3
  BC+1
}
BC0.cls <- renumerate(BC0.cls)
BCBig.cls <- renumerate(BCBig.cls)


CLS <- BC0.cls
DCOR <- BC0.dcor

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


### expand CLS matrix into glyphs
0, 1,2,3, 4, 5,6,7

to.glyph <- function(c) {
  r <- NaN
  if(c == 0)  # empty
    r <- matrix(c(0,0,0,0), nrow=2)
  if(c == 1)  # hih
    r <- matrix(c(1,1,1,0), nrow=2)
  if(c == 2)  # pc
    r <- matrix(c(0,2,2,0), nrow=2)
  if(c == 3)  # lil
    r <- matrix(c(0,3,3,3), nrow=2)
  if(c == 4)  # un
    r <- matrix(c(4,4,4,4), nrow=2)
  if(c == 5)  # hil
    r <- matrix(c(5,5,0,5), nrow=2)
  if(c == 6)  # nc
    r <- matrix(c(6,0,0,6), nrow=2)
  if(c == 7)  # lih
    r <- matrix(c(7,0,7,7), nrow=2)
  r
}
expand.cls <- function(CLS) {
  G <- mat.or.vec(nrow(CLS)*2, ncol(CLS)*2)
  for(i in 0:(nrow(CLS)-1))
    for(j in 0:(ncol(CLS)-1))
      G[(i*2+1):(i*2+2),(j*2+1):(j*2+2)] <- to.glyph(CLS[i+1,j+1])
  G
}
