# running HGC in existing Seurat object, and save the tree in the object
FindClusteringTree <- function(object, graph.type = "SNN"){
    require(Seurat)
    require(Matrix)

    SeuratDefaultAssay <- Seurat::DefaultAssay(object)

    if(graph.type == "SNN"){
        SNN_name = paste(SeuratDefaultAssay, "_snn", sep = "")
        G = Seurat::as.sparse(object@graphs[[SNN_name]])
        ClusteringTree = HGC.dendrogram(G)
        object@graphs[["ClusteringTree"]] = ClusteringTree
        return(object)
    }else if(graph.type == "KNN"){
        KNN_name = paste(SeuratDefaultAssay, "_nn", sep = "")
        G = Seurat::as.sparse(object@graphs[[KNN_name]])
        G = G + t(G)
        ClusteringTree = HGC.dendrogram(G)
        object@graphs[["ClusteringTree"]] = ClusteringTree
        return(object)
    }else{
        stop("The kind of graph is not used by Seurat.")
    }
}

# running HGC in existing igraph object, and output the clustering tree
# we will try to save the tree in the igraph object in later version
cluster_HGC <- function(graph){
    require(igraph)
    G <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
    G.ClusteringTree <- HGC.dendrogram(G = G)
    return(G.ClusteringTree)
}