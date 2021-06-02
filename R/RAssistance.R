# The function is similar to the dendrogram in "scipy.cluster.hierarchy".
# The function output the order of leaves in the clustering tree.

# The function is not open to users.

get.LeavesOrder <- function(linkageMatrix){
    return(get_leaves(linkageMatrix))
}

