# The function is similar to the dendrogram in "scipy.cluster.hierarchy".
# The function output the order of leaves in the clustering tree.

# The function is not open to users.

get.LeavesOrder <- function(linkageMatrix){
    return(get_leaves(linkageMatrix))
}

# Check the input G as dgCMatrix
# The function is not open to users.

checkGraph <- function(G){
    if(class(G)[1] != "dgCMatrix"){
        stop("Wrong graph data structure")
    }
    if(nrow(G) != ncol(G)){
        stop("Wrong graph data structure")
    }
    if(min(G) < 0){
        stop("Wrong graph data structure")
    }
    if(sum(is.na(G)) > 0){
        stop("Wrong graph data structure")
    }
    if(sum(is.infinite(G)) > 0){
        stop("Wrong graph data structure")
    }
}

# Check the input integer parameter
# The function is not open to users.

checkInteger <- function(t.parameter, str = "k"){
    if(t.parameter < 1 || round(t.parameter) != t.parameter){
        order <- paste(str, " should be positive integer.",
                        sep = "")
        stop(order)
    }
}

# check the input data structure
# The function is not open to users.

checkDataStructure <- function(object,
                                object.name,
                                expected.class = "matrix"){
    if(class(object)[1] != expected.class){
        order <- paste("Wrong input data structure for ",
                        object.name, ".", sep = "")
        stop(order)
    }
}
