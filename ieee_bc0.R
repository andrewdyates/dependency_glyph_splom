load("BC.RData")
load("data/methyl_mrna.RData")
CLS <- as.matrix(BC0.cls)
DCOR <- as.matrix(BC0.dcor)
source("lib.R")

pdf("bc0.dcor.pdf", width=20, height=20)
heatmap.3(DCOR)
dev.off()

pdf("bc0.gsplom.nogrid.pdf", width=10, height=10)
par(mar = c(0, 0, 0, 0))
R <- splom(CLS, DCOR, asGlyphs=T, draw.labs=F, grid=F)
dev.off()

pdf("bc0.gsplom.grid.pdf", width=10, height=10)
par(mar = c(0, 0, 0, 0))
R <- splom(CLS, DCOR, asGlyphs=T, draw.labs=F, grid=T)
dev.off()

meth.ids <- rownames(CLS)
mrna.ids <- colnames(CLS)
meth.ids <- meth.ids[order.dendrogram(R$Rhclust)]
mrna.ids <- mrna.ids[order.dendrogram(R$Chclust)]
  
qq.row <- match(meth.ids, rownames(methyl))
qq.col <- match(mrna.ids, rownames(mRNA))

A <- mRNA[qq.col[20:30],]
B <- methyl[qq.row[30:40],]
C <- rbind(A,B)


pdf("bc0.splom.pdf", width=30, height=30)
pairs(t(C))
dev.off()
