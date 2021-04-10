\name{MST.Construction}
\alias{MST.Construction}
\title{
Building Unweighted Minimum Spanning Tree Graph
}
\description{
This function builds an Unweighted Minimum Spanning Tree (MST) graph 
in the input feature space using Euclidean distance metric.
}
\usage{
MST.Construction(mat)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
}
\details{
This function builds a MST graph of the input data. The main
function come from the R package \code{ape}.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input samples.
The matrix C is the adjacency matrix of the built MST graph. C[i,j] = 1  
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- MST.Construction(Pollen.PCs)
}
