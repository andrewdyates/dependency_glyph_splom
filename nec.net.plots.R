## SAMPLE USE:
##
# z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
# zz83 <- split(names(z83),z83)
# syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
# clust.coords.83 <- clust.names.to.idx(zz83,syms)
# pdf("whatever.pdf")
# G <- splom()....
# plot.rects.coords(clust.coords.83)
# map.plot.rects(clust.coords.83, SOME.HIT.MATRIX)
# dev.off()

plot.rects <- function(zz, syms) {
  n <- length(syms)
  for(i in 1:length(zz)) {
    print(paste("cluster",i,"size",length(zz[[i]])))
    select.i <- which(syms %in% zz[[i]])
    x0 <- min(select.i)
    x1 <- max(select.i)
    y0 <- n-max(select.i)+1
    y1 <- n-min(select.i)+1
    rect(x0,y0,x1,y1, col=rgb(0,0,0,0.4))
  }
}
clust.names.to.idx <- function(zz,syms) {
  R <- list()
  n <- length(syms)
  for(i in 1:length(zz)) {
    R[[i]] <- list()
    select.i <- which(syms %in% zz[[i]])
    R[[i]]$x0 <- min(select.i)
    R[[i]]$x1 <- max(select.i)
    R[[i]]$y0 <- n-max(select.i)+1
    R[[i]]$y1 <- n-min(select.i)+1
  }
  R
}
plot.rects.coords <- function(coords) {
  for(c in coords) {
    rect(c$x0,c$y0,c$x1,c$y1, col=rgb(0,0,0,0.4))
  }
}
map.plot.rects <- function(coords,map) {
  ZZZ <- which(map, arr.ind=TRUE)
  n <- dim(ZZZ)[1]
  if (n==0) return()
  for(i in 1:dim(ZZZ)[1]) {
    rect(coords[[ZZZ[i,1]]]$x0, coords[[ZZZ[i,2]]]$y0, coords[[ZZZ[i,1]]]$x1, coords[[ZZZ[i,2]]]$y1, col=rgb(0.5,0,0,0.3))
    rect(coords[[ZZZ[i,2]]]$x0, coords[[ZZZ[i,1]]]$y0, coords[[ZZZ[i,2]]]$x1, coords[[ZZZ[i,1]]]$y1, col=rgb(0.5,0,0,0.3))
  }
}

# Greedily remove rows with NA classes
remove.na.features <- function(BOOL.TF) {
  R <- list()
  R$na.counts <- c()
  R$worst.counts <- c()
  R$worst.na.rm <- c()
  # Remove features that are not PC class with itself
  pc.qq <- diag(BOOL.TF) == 2
  R$BOOL.TF.F <- BOOL.TF[pc.qq,pc.qq]
  R$not.pc.rm <- rownames(BOOL.TF)[!pc.qq]
  # while NA (0) boolean relations exist, remove feature with most NA relations
  while(sum(R$BOOL.TF.F==0)>0) {
    R$na.counts <- c(R$na.counts, sum(R$BOOL.TF.F==0))
    sum.zeros <- apply(R$BOOL.TF.F==0,1,sum)
    worst <- which.max(sum.zeros)
    qq <- 1:dim(R$BOOL.TF.F)[1]!=worst
    R$BOOL.TF.F <- R$BOOL.TF.F[qq,qq]
    R$worst.na.rm <- c(R$worst.na.rm, rownames(R$BOOL.TF.F)[worst])
    R$worst.counts <- c(R$worst.counts, sum.zeros[worst])
  }
  R
}
