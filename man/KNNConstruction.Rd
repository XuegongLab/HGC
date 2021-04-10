\name{KNN.Construction}
\alias{KNN.Construction}
\title{
Building Unweighted K Nearest Neighbor Graph
}
\description{
This function builds an Unweighted K Nearest Neighbor (KNN) graph 
in the input feature space using Euclidean distance metric.
}
\usage{
KNN.Construction(mat, k)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
    \item{k}{
    the number of nearest neighbors for building the KNN graph.
    }
}
\details{
This function builds a KNN graph of the input data. The main function 
comes from the R package \code{RANN}.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input 
samples.
The matrix C is the adjacency matrix of the built KNN graph. C[i,j] = 1  
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- KNN.Construction(Pollen.PCs)
}
