\name{RNN.Construction}
\alias{RNN.Construction}
\title{
Building Unweighted \eqn{\epsilon} Nearest Neighbor Graph
}
\description{
This function builds an \eqn{\epsilon} Nearest Neighbor graph in the
input feature space using Euclidean distance metric.
}
\usage{
RNN.Construction(mat, max_dist)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
    \item{max_dist}{
    the threshold distance. The edges whose lengths are less than
    \code{max_dist} will be kept in the graph.
    }
}
\details{
The function builds an \eqn{\epsilon} Nearest Neighbor graph  which
saved as a sparse matrix.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input
samples.
The matrix C is the adjacency matrix of the built RNN graph. C[i,j] = 1
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- RNN.Construction(Pollen.PCs, 20)
}
