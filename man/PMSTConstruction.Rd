\name{PMST.Construction}
\alias{PMST.Construction}
\title{
Building Unweighted Perturbed Minimum Spanning Tree Graph
}
\description{
This function builds an Unweighted Perturbed Minimum Spanning Tree
(PMST) graph in the input feature space using Euclidean distance metric.
}
\usage{
PMST.Construction(mat, iter, r)
}
\arguments{
    \item{mat}{
    the input data saved as a numerical matrix. The columns are the
    features and the rows are the samples.
    }
    \item{iter}{
    the number of perturbation.
    }
    \item{r}{
    the parameter about the strength of the perturbation.
    }
}
\details{
The function builds a PMST graph of the input data. PMST is the
combination of a number of MSTs, which are built in the pertured data
spaces.
}
\value{
An n by n binary dgCMatrix object C, where n is the number of input
samples.
The matrix C is the adjacency matrix of the built PMST graph. C[i,j] = 1
means that there is an edge between sample i and sample j.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- PMST.Construction(Pollen.PCs)
}
