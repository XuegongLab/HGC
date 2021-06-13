\name{FindClusteringTree}
\alias{FindClusteringTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The HGC algorithm embedded in Seurat pipeline
}
\description{
The function runs hierarchical clustering with \code{HGC.dendrogram} 
on the SNN or KNN calculated by the \code{Seurat} 
pipeline. The output clustering tree is also packaged in the 
\code{Seurat} object.
}
\usage{
FindClusteringTree(object, graph.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    The Seurat object containing the graphs built with scRNA-seq data.
}
    \item{graph.type}{
    The type of graphs used for the hierarchical clustering, could be 
    "SNN" or "KNN". The default value is "SNN".
}
}
\details{
For the KNN graph, we symmetrize it by adding its transposition on 
the graph. And for the details of data preprocessing and graph construction 
by \code{Seurat}, please check the \code{Seurat} vignettes.
}
\value{
An \code{Seurat} object. The clustering tree is saved under the item 
\code{graphs}, i.e. \code{object@graphs$ClusteringTree}.
}

\note{
The function needs the R package \code{Seurat}. We recommend that the 
version of \code{Seurat} is higher than version 3.0.
}

\examples{
## Do not run
# require(Seurat)

## DemoData is a input gene expression matrix.
# DemoData.seuratobj <- CreateSeuratObject(counts = DemoData,
#                                          min.cells = 20)
# DemoData.seuratobj <- NormalizeData(object = DemoData.seuratobj,
#                                     verbose = F)
# DemoData.seuratobj <- ScaleData(object = DemoData.seuratobj,
#                                 features = row.names(DemoData.seuratobj),
#                                 verbose = F)
# DemoData.seuratobj <- FindVariableFeatures(object = DemoData.seuratobj,
#                                            nfeatures = 2000, verbose = F)
# DemoData.seuratobj <- RunPCA(object = DemoData.seuratobj,
#                              npcs = 100, verbose = F)
# DemoData.seuratobj <- FindNeighbors(object = DemoData.seuratobj,
#                                     nn.eps = 0.5, k.param = 30,
#                                     dims = 1:25, verbose = F)
# DemoData.seuratobj <- FindClusteringTree(object = DemoData.seuratobj,
#                                          graph.type = "SNN")
}
