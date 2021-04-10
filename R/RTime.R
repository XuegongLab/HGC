# HGC.time is a function to observe the nearest neighbor chain length and the
# average neighbor number during the clustering.

HGC.parameter <- function(G){
    # Run clustering
    checkGraph(G)
    record = HierarCluster_paris_time(G)
    rownames(record) = c("ChainLength", "AvgNeighNum")

    return(record)
}
