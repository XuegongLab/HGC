# HGC.time is a function to observe the nearest neighbor chain length and the
# average neighbor number during the clustering.

HGC.time <- function(G){
  # Run clustering
  record = HierarCluster_paris_time(G)
  rownames(record) = c("ChainLength", "AvgNeighNum")

  return(record)
}