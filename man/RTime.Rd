\name{HGC.time}
\alias{HGC.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recording the Parameters of the Graph-based Hierarchical Clustering
}
\description{
The length of the nearest neighbor chain and the average neighbor number are related to the running time of the \code{HGC.paris} algorithm. The function records and outputs the values in each iteration to infer the time complexity of \code{HGC.paris}.
}
\usage{
HGC.time(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
  an undirected graph saved as a dgCMatrix. The matrix G is the adjacency matrix of the graph, and element G[i,j] is the weight of the edge between node i and node j. Zeros in the matrix mean no link between nodes here.
}
}
\details{
The function finishes the whole clustering progress as \code{HGC.paris}. So it may spend same time.
}
\value{
A 2 by m matrix. The two rows of the matrix are the NN chain length and the average neighbor number. m is equal to n-s, where s is the number of unconnected parts in the graph.
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
