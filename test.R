install.packages("HGC_0.1.0.tar.gz", repos = NULL, type = "source")
library(HGC)
library(Matrix)

G = matrix(1:25,5,5)
G = G + t(G)
G = as(G, "dgCMatrix")
tree = HGC.paris(G)
tree

record = HGC.time(G)
record