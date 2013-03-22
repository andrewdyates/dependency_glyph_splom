# COMPUTE DISTANCE MATRIX
# ====================
# Define glyph pairwise distance (from original UNL==0 enumeration)

#CLS.ENUM <- list("0"="NA", "1"="HIH", "2"="PC", "3"="LIL", "4"="UNL", "5"="HIL", "6"="NC", "7"="LIH", "8"="PAD")
# 0 = NA
# 0:na, 1:hih, 2:pc, 3:lil, 4:unl, 5:hil, 6:nc, 7:lih

dist.glyph <- matrix(c(r1,r2,r3,r4,r5,r6,r7,r8), 8,8)




# Compute distance between two vectors of glyphs
glyph.f <- function(A,B) {
  # glyphs are enumerated from zero already
  sqrt(sum(dist.glyph[as.numeric(A*8+(B+1))]**2))
}

# Compute distance matrix (this is rather slow, but faster than other options
gen.glyph.dist.m <- function(BC) {
  D <- matrix(0,nrow(BC),nrow(BC))
  for (i in 1:nrow(BC)) {
    for (j in (i+1):nrow(BC)) {
      if (j <= nrow(BC)) {
        D[j,i] <- glyph.f(BC[i,], BC[j,])
      }
    }
  }
  as.dist(D)
}




