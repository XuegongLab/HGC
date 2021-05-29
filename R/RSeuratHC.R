FindClusteringTree <- function(object, graph.type = "SNN", random.seed = 0){
  set.seed(random.seed)
  require(Seurat)
  require(Matrix)

  SeuratDefaultAssay <- DefaultAssay(object)

  if(graph.type == "SNN"){
    SNN_name = paste(SeuratDefaultAssay, "_snn", sep = "")
    G = as.sparse(object@graphs[[SNN_name]])
    ClusteringTree = HGC.paris(G)
    object@graphs[["ClusteringTree"]] = ClusteringTree
    return(object)
  }else if(graph.type == "KNN"){
    KNN_name = paste(SeuratDefaultAssay, "_nn", sep = "")
    G = as.sparse(object@graphs[[KNN_name]])
    G = G + t(G)
    ClusteringTree = HGC.paris(G)
    object@graphs[["ClusteringTree"]] = ClusteringTree
    return(object)
  }else{
    stop("The kind of graph is not used by Seurat.")
  }

}
