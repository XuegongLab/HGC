# HGC: A graph-based hierarchical clustering algorithm


## Introduction

**`HGC`** is an R package for **HGC: A graph-based hierarchical clustering method for scRNA-seq data**. It could provide hierarchical clustering result for large scale single cell RNA sequencing (scRNA-seq) datasets.

**`HGC`** catch the hierarchical information in the undirected weighted graphs. For scRNA-seq data, the clustering result is based on the shared nearest neighbor (SNN) graphs built in feature spaces. For scRNA-seq data, we organize a pipeline around HGC, and its overview is shown here.

![overview](https://github.com/stevenhuakui/HGC/blob/master/figures/overview.png)

[//]: #  (## Citation)

[//]: #  (If you use **`HGC`** in published research, please cite:)


## Installation

To install the *developmental version* from [**GitHub**](https://github.com/stevenhuakui/HGC/):

```{r Installation from GitHub, eval = FALSE}
if(!require(devtools)) install.packages("devtools")
devtools::install_github("stevenhuakui/HGC", build_vignettes = TRUE)
```

To load the installed **`HGC`** in R:

```{r Load HGC, eval = FALSE}
library(HGC)
```

The `HGC` package needs the support of CPP11. Please check your system environment.


## Input

**`HGC`** takes one inputs: the graph `G`.

The input `G` is the adjacency matrix of the graph. It should be stored as the dgCMatrix supported by the R package `Matrix`. The row names and column names of the matrix are both the nodes' names, that are the cells' names for scRNA-seq data. The element of `G` are the weights of edges in the graph. Zeros in `G` mean no link between the correspoding nodes.


## Usage

### With dgCMatrix input

Here is an example to run **`HGC`** with random input:

```{r demo1, eval = FALSE}
require(Matrix)

G = matrix(1:25,5,5)
G = G + t(G)
G = as(G, "dgCMatrix")
tree = HGC.paris(G)
record = HGC.time(G)
```

### With inputs from different scRNA-seq pipelines

There are many popular scRNA-seq pipelines building and utilizing graphs during the analysis. We try to collect and organize a guide to apply HGC in the exisiting graphs in the pipelines. The guide is stored in [**GitHub**](https://github.com/stevenhuakui/HGC_plot), as "Loading_Different_Graphs.ipynb".


## Output

The clustering tree by HGC is stored as `hclust` object supported by the R package `stats`. The detailed description is in the help pages of `HGC.paris` or `hclust`.


## Visualization of results

We collect many exisiting methods to visualize the clustering trees in R. There are some examples.

<table><tr>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig1.png" width="256"/>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig2.png" width="256"/>
    <td><img src="https://github.com/stevenhuakui/HGC/blob/master/figures/fig3.png" width="256"/>
</tr></table>

The detailed codes for plotting are stored in [**GitHub**](https://github.com/stevenhuakui/HGC_plot), as "Plotting_Dendrogram_And_Low_Dimension_Visualization.ipynb".


## Help

Use the following code in R to get access to the help documentation for **`HGC`**:

```{r help1, eval = FALSE}
# Documentation for HGC
?HGC.paris
```

```{r help2, eval = FALSE}
# Documentation for HGC
?HGC.time
```

You are also welcome to view and post questions in github or contact the author by email for help.