FindClusteringTree <- function(object, graph.type = "SNN", random.seed = 0){
  set.seed(random.seed)
  require(Seurat)
  require(Matrix)

  if(graph.type == "SNN"){
    G = as.sparse(object@graphs$RNA_snn)
    ClusteringTree = HGC.paris(G)
    object@graphs[["ClusteringTree"]] = ClusteringTree
    return(object)
  }else if(graph.type == "KNN"){
    G = as.sparse(object@graphs$RNA_nn)
    G = G + t(G)
    ClusteringTree = HGC.paris(G)
    object@graphs[["ClusteringTree"]] = ClusteringTree
    return(object)
  }else{
    stop("The kind of graph is not used by Seurat.")
  }

}
