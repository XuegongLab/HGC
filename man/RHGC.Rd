\name{HGC.paris}
\alias{HGC.paris}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graph-based Hierarchical Clustering
}
\description{
Hierarchical clustering on a given undirected graph.
}
\usage{
HGC.paris(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
  an undirected graph saved as a dgCMatrix. The matrix G is the adjacency matrix of the graph, and element G[i,j] is the weight of the edge between node i and node j. Zeros in the matrix mean no link between nodes here.
}
}
\details{
The function runs a hierarchical clustering on the given graph. The node pair sampling ratio is calculated among the connected nodes in the graph as the distance. And the algorithm will merge the closest nodes together iteratively, untill there is only one node left in the graph.
}
\value{
An object of class \code{hclust} defined by the \code{hclust} function in the \code{stats} package. It is a list containing the clustering tree information with the components:
\item{merge }{an \emph{n-1} by 2 matrix. It records the two nodes in each merging step.}
\item{height }{a set of \emph{n-1} real values. It is the height of the non-leaf nodes in the tree.}
\item{order }{a vector giving the permutation of the original observations suitable for plotting.}
\item{labels }{labels for the objects being clustered. Same as the rownames of G in default.}
\item{call }{the call which produced the result.}
\item{method }{the cluster method that has been used.}
\item{dist.method }{the distance used here.}

More details about the components are in the \link[stats]{hclust}.
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(Matrix)

G = matrix(1:25,5,5)
G = G + t(G)
G = as(G, "dgCMatrix")
tree = HGC.paris(G)
record = HGC.time(G)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
