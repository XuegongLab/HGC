# HGC: fast hierarchical clustering for large-scale single-cell data 


## Introduction

**`HGC`** is an R package for conducting fast hierarchical clustering on large-scale single-cell RNA sequencing (scRNA-seq) datasets.

The basic idea of **`HGC`** is to perform hierarchical clustering on the shared nearest neighbor (SNN) graphs of cells. The overview of HGC is shown as follows.

![overview](https://github.com/stevenhuakui/HGC/blob/master/figures/overview2.png)

[//]: #  (## Citation)

[//]: #  (If you use **`HGC`** in published research, please cite:)


## Installation

To install the *developmental version* from [**GitHub**](https://github.com/stevenhuakui/HGC/):

```{r Installation from GitHub, eval = FALSE}
if(!require(devtools))
  install.packages("devtools")
devtools::install_github("stevenhuakui/HGC", build_vignettes = TRUE)
```

To load the installed **`HGC`** in R:

```{r Load HGC, eval = FALSE}
library(HGC)
```

The `HGC` package needs the support of CPP11. Please check your system environment.





## Usage

The core function of **`HGC`** takes a graph `G` as the input. Here `G` is the adjacency matrix of a graph of cells, where the element in `G` are the weights of edges in the graph, and zero means no link between the corresponding node pairs. `G` should be a dgCMatrix supported by the R package `Matrix`.

**`HGC`** could be seamlessly used in the popular `Seurat` pipeline. It also accepts self-build graphs or cell graph from other scRNA-seq pipelines. Here are some examples of the usage.

### With Seurat object as input

For analysis in `Seurat` pipeline, use the `FindClusteringTree` function:

```{r demo2, eval = FALSE}
# Use DemoData to represent an gene expression matrix
require(Seurat)

DemoData.seuratobj <- CreateSeuratObject(counts = DemoData,
                                         min.cells = 20)
DemoData.seuratobj <- NormalizeData(object = DemoData.seuratobj,
                                    verbose = F)
DemoData.seuratobj <- ScaleData(object = DemoData.seuratobj,
                                features = row.names(DemoData.seuratobj),
                                verbose = F)
DemoData.seuratobj <- FindVariableFeatures(object = DemoData.seuratobj,
                                           nfeatures = 2000, verbose = F)
DemoData.seuratobj <- RunPCA(object = DemoData.seuratobj,
                             npcs = 100, verbose = F)
DemoData.seuratobj <- FindNeighbors(object = DemoData.seuratobj,
                                    nn.eps = 0.5, k.param = 30,
                                    dims = 1:25, verbose = F)
DemoData.seuratobj <- FindClusteringTree(object = DemoData.seuratobj,
                                         graph.type = "SNN")
```

### With dgCMatrix as input

For other graphs, use the function `HGC.paris`:

```{r demo1, eval = FALSE}
require(Matrix)

G = matrix(1:25,5,5)
G = G + t(G)
G = as(G, "dgCMatrix")
tree = HGC.paris(G)
record = HGC.time(G)
```



### With graphs from other scRNA-seq pipelines

There are many popular scRNA-seq pipelines building and utilizing graphs in their analysis. We try to collect and organize a guide to apply **`HGC`** in the exisiting graphs in the pipelines. The guide is stored in [**GitHub**](https://github.com/stevenhuakui/HGC/tree/master/HGC_plotting_guide), as "Loading_Different_Graphs.ipynb".


## Output

The dendrogram given by **`HGC`** is stored as `hclust` object supported by the R package `stats`. The detailed descriptions can be found in the help pages of `HGC.paris` or `hclust`.


## Visualization of results

We collect many exisiting methods to visualize the clustering trees in R. Here are some examples.

<table><tr>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig1.png" width="256"/>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig2.png" width="256"/>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig3.png" width="256"/>
</tr></table>

The detailed codes for plotting are stored in [**GitHub**](https://github.com/stevenhuakui/HGC/tree/master/HGC_plotting_guide), as "Plotting_Dendrogram_And_Low_Dimension_Visualization.ipynb".


## Help

Use the following code to get access to the help documentation of **`HGC`**:

```{r help3, eval = FALSE}
# Documentation for HGC
?FindClusteringTree
```

```{r help1, eval = FALSE}
# Documentation for HGC
?HGC.paris
```

```{r help2, eval = FALSE}
# Documentation for HGC
?HGC.time
```

You are welcome to view and post questions in github or contact the author by email for help.
