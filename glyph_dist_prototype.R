# COMPUTE DISTANCE MATRIX
# ====================
# Define glyph pairwise distance (from original UNL==0 enumeration)

#CLS.ENUM <- list("0"="NA", "1"="HIH", "2"="PC", "3"="LIL", "4"="UNL", "5"="HIL", "6"="NC", "7"="LIH", "8"="PAD")
# 0 = NA
# 0:na, 1:hih, 2:pc, 3:lil, 4:unl, 5:hil, 6:nc, 7:lih

load("mar21.celegans.gold.dcor.cls.wtms.RData")

G.DIST.TABLE <- read.table("glyph.dists.csv", as.is=T, sep=",")

glyph.dist.f <- function(A,B) {
  # Requires that NA is mapped to 8 rather than 0
  sum(apply(cbind(A,B), 1, function(p) G.DIST.TABLE[p[1],p[2]]))
}

A <- c(1,5,7,3,5,6) # = from A: 0
B <- c(2,4,3,0,0,2) # = from A: 12
C <- c(2,4,3,0,5,2) # = from A: 11
D <- c(6,5,7,3,5,6) # = from A: 3
M <- rbind(A,B,C,D)
#      1 1 2 2 2 4 = 12
# Compute distance between two vectors of glyph enumerations.
# CLS.ENUM <- list("0"="NA", "1"="HIH", "2"="PC", "3"="LIL", "4"="UNL", "5"="HIL", "6"="NC", "7"="LIH", "8"="PAD")


# Compute distance matrix (this is rather slow, but faster than other options
gen.glyph.dist.m <- function(BC, recast.na.0=T) {
  if(recast.na.0)
    BC[BC==0] <- 8
  n <- nrow(BC)
  as.dist(outer(1:n,1:n, FUN = Vectorize( function(i,j) glyph.dist.f(BC[i,],BC[j,]) )))
}

