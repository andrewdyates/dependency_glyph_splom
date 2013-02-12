## Generate glyphs without coloring
## Custom Color Scales
source("lib.R")
load("BC.RData") ## generated from data/*.tab files and compiled using code in prototype.R
CLS <- BC0.cls
DCOR <- BC0.dcor

asGlyphs=FALSE; pad=FALSE; N=15; MIN=0.1; MAX=0.8; MOST=1; LEAST=0; DCOR.weight=2; useRaster=FALSE


