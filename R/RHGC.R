# HGC.paris is the main function to cluster. It outputs the clustering tree.

HGC.dendrogram <- function(G){
    # Run clustering
    checkGraph(G)
    hg.linkageMatrix = HierarCluster_paris(G)

    # check the linkageMatrix
    hg.linkageMatrix2 <- hg.linkageMatrix
    hg.linkageMatrix.temp <- as.vector(hg.linkageMatrix[,3])

    # The unconnected branches problem
    link.na <- which(is.infinite(hg.linkageMatrix.temp))
    if(length(link.na) > 0){
        hg.linkageMatrix.temp <- hg.linkageMatrix.temp[-link.na]
        hg.linkageMatrix2[which(is.infinite(hg.linkageMatrix))] <- 10*max(
            hg.linkageMatrix.temp, na.rm = TRUE)
    }

    hg.tree2 = list()
    class(hg.tree2) <- 'hclust'
    hg.tree2$call = NA
    hg.tree2$dist.method = "Node pair sampling ratio"
    hg.tree2$method = "HGC"

    hg.tree2$height = hg.linkageMatrix2[,3]

    N = ncol(G)
    merge.m = hg.linkageMatrix2[,seq_len(2)] + 1
    merge.m[merge.m<=N] = merge.m[merge.m<=N]*-1
    merge.m[merge.m>N] = merge.m[merge.m>N]-N

    hg.tree2$merge = merge.m

    hg.tree2$order = get.LeavesOrder(hg.linkageMatrix2) + 1

    if(length(rownames(G)) == nrow(G)){
        hg.tree2$labels = rownames(G)
    }else{
        hg.tree2$labels = NA
    }

    return(hg.tree2)
}
