\name{SNN.Construction}
\alias{SNN.Construction}
\title{
Building Unweighted Shared Nearest Neighbor Graph
}
\description{
This function builds a Shared Nearest Neighbor (SNN) graph in the input 
feature space using Euclidean distance metric.
}
\usage{
SNN.Construction(mat, k, threshold)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
    \item{k}{
    the number of nearest neighbor number to build the original KNN.
    }
    \item{threshold}{
    the threshold parameter for the Jaccard index. The edges in KNN
    whose Jaccard indices are lower than it will be removed in building 
    the SNN.
    }
}
\details{
The function builds an SNN which saved as a sparse matrix.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input 
samples.
The matrix C is the adjacency matrix of the built SNN graph. C[i,j] = 1  
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- SNN.Construction(Pollen.PCs)
}
