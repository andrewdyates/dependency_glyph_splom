source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
# loaded BC0.cls, BC0.dcor, BCBig.cls, BCBig.dcor

# Test heatmap
heatmap.3(BC0.dcor)
heatmap.3(BCBig.dcor)

# Test class clustering and plotting
CLS <- order.cls(BC0.cls)
draw.glyphs(CLS)
draw.glyphs(CLS, grid=1)

CLS.big <- order.cls(BCBig.cls)
draw.glyphs(CLS.big)

# Test glyph expansion and plotting
GLY <- expand.cls(CLS)
draw.glyphs(GLY)
draw.glyphs(GLY, grid=2)
draw.glyphs(GLY, grid=2, grid.col="#ff0000")
# Test pixel padding
PAD <- expand.cls(CLS, pad=TRUE)
draw.glyphs(PAD)
draw.glyphs(PAD, grid=3)

# Test pixel plotting convenience function.
# single pixel plot



plot.pix(CLS.big, scale=1, fname="cls.big")

plot.pix(CLS, scale=10, fname="cls.test.scale10")
plot.pix(CLS, scale=10, fname="cls.test.scale10.grid1", grid=1)
plot.pix(CLS, scale=1, fname="cls.test.scale1.pdf")
plot.pix(CLS, scale=1, fname="cls.test.scale1.grid1.png", grid=1)
# glyph plot
plot.pix(GLY, scale=10, fname="gly.test.scale10")
plot.pix(GLY, scale=10, fname="gly.test.scale10.grid2", grid=2)
plot.pix(GLY, scale=1, fname="gly.test.scale1")
plot.pix(GLY, scale=1, fname="gly.test.scale1.grid2", grid=2)
# glyph plot with padding
plot.pix(PAD, scale=10, fname="pad.test.scale10")
plot.pix(PAD, scale=10, fname="pad.test.scale10.grid3", grid=3, grid.col="#e5e5ff")
plot.pix(PAD, scale=1, fname="pad.test.scale1.")
plot.pix(PAD, scale=1, fname="pad.test.scale1.grid3", grid=3, grid.col="#e5e5ff")

# Test vector plotting convenience function.
plot.vtr(CLS, fname="cls.test")
plot.vtr(CLS, fname="cls.test.grid1", grid=1)
plot.vtr(CLS, width=20, height=20, fname="cls.test.grid1.w20.h20", grid=1)
plot.vtr(CLS, width=20, fname="cls.test.w20.pdf")
plot.vtr(CLS, width=20, fname="cls.test.w20.grid1.png", grid=1)
plot.vtr(CLS, width=3, fname="cls.test.w3.grid1.png", grid=1)
# glyph plot
plot.vtr(GLY, fname="gly.test")
plot.vtr(GLY, fname="gly.test.grid2", grid=2)
plot.vtr(GLY, fname="gly.test.grid2.ltblue", grid=2, grid.col="#e5e5ff")
plot.vtr(GLY, width=20, height=20, fname="gly.test.grid2.w20.h20", grid=2)
plot.vtr(GLY, width=3, fname="gly.test.grid2.w3", grid=2)
# pad plot
plot.vtr(PAD, fname="pad.test")
