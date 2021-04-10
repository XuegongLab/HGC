\name{HGC.parameter}
\alias{HGC.parameter}
\title{
Recording the Parameters of the Graph-based Hierarchical Clustering
}
\description{
This function records and outputs the length of the nearest neighbor 
chain and the average neighbor
number in each iteration of hierarchical clustering. These values can 
be used for the time complexity analysis of \code{HGC.dendrogram}.
}
\usage{
HGC.parameter(G)
}
\arguments{
    \item{G}{
    an undirected graph saved as a dgCMatrix. The matrix G is the adjacency
    matrix of the graph, and element G[i,j] is the weight of the edge between
    node i and node j. Zeros in the matrix mean no link between nodes here.
    }
}
\details{
This function contains the whole function of \code{HGC.dendrogram}, but will 
record the key parameters during the whole process.  The function is provided 
for advanced users to conduct time complexity analysis on their own data. 
The construction of the dendrogram is a recursive procedure of two steps: 
1.finding the nearest neighbour pair, 2. merge the node pair and update the 
graph. For different data structures of graph, there's a trade-off between 
the time consumptions of the two steps. Generally speaking, storing more 
information about the graph makes it faster to find the nearest neighbour 
pair (step 1) but slower to update the graph (step 2). We have experimented 
serval datasets and chosen the best data structure in \code{HGC.dendrogram} 
for the overall efficiency.
}
\value{
A 2 by m matrix. The two rows of the matrix are the nearest neighbor chain 
length and the average neighbor number. m is equal to n-s, where s is 
the number of unconnected parts in the graph.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
G <- SNN.Construction(Pollen.PCs, 25, 0.15)
record = HGC.parameter(G)
}
