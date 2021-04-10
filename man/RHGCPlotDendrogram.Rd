\name{HGC.PlotDendrogram}
\alias{HGC.PlotDendrogram}
\title{
Visualizing the dendrogram
}
\description{
The function will plot the dendrogram, with 
different colors for different clusters.
}
\usage{
HGC.PlotDendrogram(tree, k, plot.label, labels)
}
\arguments{
    \item{tree}{
    the input clustering tree saved as \code{hclust} data structure.
    }
    \item{k}{
    the number of clusters to cut the tree into.
    }
    \item{plot.label}{
    a bool variable. It decides whether the function will add color bars.
    }
    \item{labels}{
    a data frame or a matrix to store the label information.
    Different labels should be in different columns and the
    users should name the columns correspondingly. The label
    information will show in the figure as color
    bars below the dendrogram, and each label takes one color bar.
    }
}
\details{
The function plots the clustering tree, with alternative colors showing
the clustering results and the label information. It is based on the R
package \code{dendextend} which contains contains many parameters for the
visualization. For users' convenience, most of the parameters are set to
be the default value. Advanced users could
visit the vignette of \code{dendextend} for more flexible visualization.
}
\value{
The function will return 1 if the dendrogram is sucessfully drawn.
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
HGC.PlotDendrogram(tree = Pollen.ClusteringTree,
                    k = 5, plot.label = TRUE,
                    labels = Pollen.labels)
}
