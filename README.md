# HGC: fast hierarchical clustering for large-scale single-cell data
## Introduction

`HGC` (short for Hierarchical Graph-based Clustering) is an R package for 
conducting hierarchical clustering on large-scale single-cell RNA-seq 
(scRNA-seq) data. The key idea is to construct a dendrogram of cells on 
their shared nearest neighbor (SNN) graph. `HGC` provides functions for 
building cell graphs and for conducting hierarchical clustering on the graph. 
Experiments on benchmark datasets showed that `HGC` can reveal the 
hierarchical structure underlying the data, achieve state-of-the-art 
clustering accuracy and has better scalability to large single-cell data. 
For more information, please refer to the preprint of `HGC` on 
[bioRxiv](https://doi.org/10.1101/2021.02.07.430106).

`HGC.core` package is the skeleton of `HGC` and it constains only clustering functions.

## Installation

`HGC.core` could be installed from Github.

```{r Github install, eval = FALSE}
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("XuegongLab/HGC", ref = "HGC_core")
```

## Quick Start

### Input data and Run HGC
The hierarchical clustering function `HGC.dendrogram` needs the adjacency matrix as input which saved as `dgCMatrix`.

```{r, message=FALSE, warning=FALSE}
library(HGC.core)
library(Matrix)

G <- matrix(c(1:100), nrow = 10, ncol = 10)
G <- G + t(G)
G <- as(G, "dgCMatrix")
tree = HGC.dendrogram(G)
```