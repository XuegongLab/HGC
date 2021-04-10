\name{HGC.dendrogram}
\alias{HGC.dendrogram}
\title{
Hierarchical Graph-based Clustering
}
\description{
Hierarchical clustering on a given undirected graph.
}
\usage{
HGC.dendrogram(G)
}
\arguments{
    \item{G}{
    an undirected graph saved as a dgCMatrix. The matrix G is the
    adjacency matrix of the graph, where G[i,j] is the weight
    of the edge between node i and node j, and zero means
    no link.
    }
}
\details{
The function runs a hierarchical clustering on the given graph. It is a
recursive procedure of two steps, first, the node
pair sampling ratio is used as the distance metric to search the nearest
neighbor pairs. Then the neighbor pair are merged and the graph is updated.
The whole procedure is accelerated using the nearest neighbor chain algorithm.
The algorithm stops when there's only one node left in the updated graph.
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
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- SNN.Construction(Pollen.PCs, 25, 0.15)
tree = HGC.dendrogram(G)
}
