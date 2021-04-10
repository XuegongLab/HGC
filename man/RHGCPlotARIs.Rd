\name{HGC.PlotARIs}
\alias{HGC.PlotARIs}
\title{
Calculating and Visualizing ARIs of the clustering results with given labels
}
\description{
The function cut the dendrogram into specific clusters at different levels
and compared the clusterings with given labels using Adjusted Rand Index (ARI)
}
\usage{
HGC.PlotARIs(tree, k.min, k.max, labels, return.ARI)
}
\arguments{
    \item{tree}{
    the input clustering tree saved as \code{hclust} data structure.
    }
    \item{k.min}{
    the minimum number to cut the tree.
    }
    \item{k.max}{
    the maximum number to cut the tree.
    }
    \item{labels}{
    a data frame or a matrix to store the label information. Different
    labels should be in different columns and the users should name the
    columns correspondingly.
    }
    \item{return.ARI}{
    a bool variable to choose whether output the ARI matrix.
    }
}
\details{
ARI is a widely used index to evaluate the consistence between two partitions
of the same samples. This function will first cut a given tree into specific
number of clusters using the function \code{cutree}. Then it calculates the
ARIs between the clustering result and the given labels
with the help of R package \code{mclust}. The function does such cutting and
calculation for different ks between k.min and k.max.Finally it visualize
these results using a line chart. ARIs with different labels are shown as
different lines with different colors in the figure.
}
\value{
A line chart will be drawn and a matrix of the ARIs will be returned.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
Pollen.Label.Tissue <- Pollen[["Tissue"]]
Pollen.Label.CellLine <- Pollen[["CellLine"]]

Pollen.SNN <- SNN.Construction(Pollen.PCs)
Pollen.ClusteringTree <- HGC.dendrogram(G = Pollen.SNN)
Pollen.labels <- data.frame(Tissue = Pollen.Label.Tissue,
                            CellLine = Pollen.Label.CellLine)
HGC.PlotARIs(tree = Pollen.ClusteringTree,
                            k.min = 2, k.max = 15,
                            labels = Pollen.labels)
}
