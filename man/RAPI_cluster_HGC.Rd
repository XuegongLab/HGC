\name{cluster_HGC}
\alias{cluster_HGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The HGC algorithm embeded in igraph pipeline
}
\description{
The \code{HGC} with the graph saved in the \code{igraph} object. The 
output is the clustering tree saved as \code{hclust} object.
}
\usage{
cluster_HGC(graph = G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
  The graph to analyze.
}
}
\details{
The graph saved in igraph object should be undirected.
}
\value{
An object of class \code{hclust} defined by the \code{hclust} function
in the \code{stats} package. It is a list containing the clustering tree
information with the components:
\item{merge }{an \emph{n-1} by 2 matrix. It records the two nodes in
each merging step.}
\item{height }{a set of \emph{n-1} real values. It is the height of the
non-leaf nodes in the tree.}
\item{order }{a vector giving the permutation of the original
observations suitable for plotting.}
\item{labels }{labels for the objects being clustered. Same as the
rownames of G in default.}
\item{call }{the call which produced the result.}
\item{method }{the cluster method that has been used.}
\item{dist.method }{the distance used here.}

More details about the components are in the \link[stats]{hclust}.
}

\note{
The function needs the R package \code{igraph}.
}

\examples{
## Do not run
# require(igraph)
# g <- sample_pa(100, m = 2, directed = FALSE)
# G.ClusteringTree <- cluster_HGC(g)
}
