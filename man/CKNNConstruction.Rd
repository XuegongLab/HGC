\name{CKNN.Construction}
\alias{CKNN.Construction}
\title{
Building Unweighted Continuous K Nearest Neighbor Graph
}
\description{
This function builds a Continuous K Nearest Neighbor (CKNN) graph 
in the input feature space using Euclidean distance metric.
}
\usage{
CKNN.Construction(mat, k, delta)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
    \item{k}{
    the number of nearest neighbors for building the CKNN graph.
    }
    \item{delta}{
    the parameter related with the distance threshold.
    }
}
\details{
This function fist built a KNN graph from the input data.Then the CKNN 
graph is built from the KNN graph. For node i and node j in
the KNN graph, CKNN will link them if the distance d(i,j) between node
i and node j is less than \eqn{\delta} times of the geometric mean of
d_k(i) and d_k(j). Here \eqn{\delta} is the parameter, d_k(i) and d_k(j)
are the distances from node i or node j to their k nearest neighbor.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input samples.
The matrix C is the adjacency matrix of the built CKNN graph. C[i,j] = 1  
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- CKNN.Construction(Pollen.PCs)
}
