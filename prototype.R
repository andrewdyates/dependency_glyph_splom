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

# plot corresponding glyphs
rowInd <- order.dendrogram(Rhclust)
colInd <- order.dendrogram(Chclust)
CLS.dcor.ord <- CLS[rowInd, colInd]
draw.glyphs(CLS.dcor.ord)
draw.glyphs(expand.cls(CLS.dcor.ord), grid=2)
# ------------------------------



CLS <- as.matrix(BC0.cls)
DCOR <- as.matrix(BC0.dcor)
# ==============================
# COMBINED GLYPH SPLOM, NO SCALING
# --------------------


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

CLS.ord <- CLS[rowInd, colInd]
draw.glyphs(CLS.ord)
draw.glyphs(expand.cls(CLS.ord))
draw.glyphs(CLS.ord, grid=1)
draw.glyphs(expand.cls(CLS.ord), grid=2)

