# Some published graph construction methods

matrix2dgCMatrix <- function(m){
    numrow <- nrow(m)
    numcolumn <- ncol(m)
    new.dgC <- Matrix::sparseMatrix(i = rep(seq_len(numrow), times = numcolumn),
                                    j = rep(seq_len(numcolumn), each = numrow),
                                    x = as.vector(m), dims =
                                    c(numrow, numcolumn))
    return(new.dgC)
}

KNN.Construction <- function(mat, k = 20){
    checkInteger(k, "k")
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    if(k >= ncell){
        stop("k should be smaller than cell number.")
    }
    tmp <- RANN::nn2(mat, mat, k + 1, searchtype = "standard")
    neighborMatrix <- tmp[[1]][, -1]
    distMatrix <- tmp[[2]][, -1]
    knn.G <- Matrix::sparseMatrix(i = rep(seq_len(ncell), k),
                                j = as.vector(neighborMatrix),
                                x = 1, dims = c(ncell, ncell))
    knn.G <- (knn.G + Matrix::t(knn.G))/2
    knn.G[Matrix::which(knn.G>0)] = 1
    #knn.G <- as.matrix(G)

    return(knn.G)
}

RNN.Construction <- function(mat, max_dist){
    if(max_dist < 0){
        stop("The threshold distance should be non-negative.")
    }
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    D <- stats::dist(mat)
    G <- as.matrix(D)
    G[which(G>max_dist)] <- 0
    G[which(G>0)] <- 1

    #rnn.G <- Matrix::Matrix(G, sparse = TRUE)
    rnn.G <- matrix2dgCMatrix(G)
    return(rnn.G)
}

CKNN.Construction <- function(mat, k = 20, delta = 1){
    checkInteger(k, "k")
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    if(k >= ncell){
        stop("k should be smaller than cell number.")
    }
    if(delta <= 0){
        stop("The parameter delta should be positive.")
    }

    Dist_original <- as.matrix(stats::dist(mat))
    ## Build the knn on the original space
    tmp <- RANN::nn2(mat, mat, k + 1, searchtype = "standard")
    neighborMatrix <- tmp[[1]][, -1]
    distMatrix <- tmp[[2]][, -1]

    cknn.G <- matrix(0, nrow = ncell, ncol = ncell)
    for(i in seq_len(ncell-1)){
        for(j in seq(i+1,ncell)){
            distk_i <- distMatrix[i,k]
            distk_j <- distMatrix[j,k]
            if(Dist_original[i,j] < delta*sqrt(distk_i*distk_j)){
                cknn.G[i,j] = 1
                cknn.G[j,i] = 1
            }
        }
    }

    cknn.G <- matrix2dgCMatrix(cknn.G)
    return(cknn.G)
}

SNN.Construction <- function(mat, k = 20, threshold = 1/15){
    checkInteger(k, "k")
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    if(k >= ncell){
        stop("k should be smaller than cell number.")
    }
    if(threshold < 0 || threshold > 1){
        stop("The parameter threshold should be within zero to one.")
    }

    tmp <- RANN::nn2(mat, mat, k + 1, searchtype = "standard")
    neighborMatrix <- tmp[[1]][, -1]
    snn.G <- ComputeSNN(nn_ranked = neighborMatrix,
                        prune = threshold)

    return(snn.G)
}

MST.Construction <- function(mat){
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    D <- stats::dist(mat)
    dmst <- ape::mst(D)

    mst.G <- matrix(0, nrow = ncell, ncol = ncell)
    mst.G[which(dmst > 0)] = 1

    mst.G <- matrix2dgCMatrix(mst.G)
    return(mst.G)
}

PMST.Construction <- function(mat, iter = 20, r = 0.5){
    checkInteger(iter, "iter")
    if(r <= 0){
        stop("The parameter r should be positive.")
    }
    ncell <- nrow(mat)
    nfeature <- ncol(mat)
    D <- stats::dist(mat)
    mst.G <- matrix(0, nrow = ncell, ncol = ncell)
    ## Build the knn on the original space
    k <- 3
    tmp <- RANN::nn2(mat, mat, k + 1, searchtype = "standard")
    neighborMatrix <- tmp[[1]][, -1]
    distMatrix <- tmp[[2]][, -1]
    # for each node, record the coordinate of the
    # average of its three nearest neighbors
    neighborLocation <- mat
    for(i in seq_len(ncell)){
        neighbor1 <- neighborMatrix[i,1]
        neighbor2 <- neighborMatrix[i,2]
        neighbor3 <- neighborMatrix[i,3]
        neighbor1.site <- as.vector(mat[neighbor1,])
        neighbor2.site <- as.vector(mat[neighbor2,])
        neighbor3.site <- as.vector(mat[neighbor3,])
        new.neighbor <- (neighbor1.site +
            neighbor2.site + neighbor3.site)/3
        neighborLocation[i,] <- new.neighbor
    }
    ## Perturbation
    for(i in seq_len(iter)){
        new.mat1 <- mat
        new.mat2 <- neighborLocation

        random.r1 <- stats::runif(n = ncell, min = 0, max = r)
        random.r2 <- rep(1, ncell) - random.r1

        new.mat <- new.mat1*random.r1 + new.mat2*random.r2
        new.D <- stats::dist(new.mat)
        new.mst <- ape::mst(new.D)

        new.mst.G <- matrix(0, nrow = ncell, ncol = ncell)
        new.mst.G[which(new.mst > 0)] = 1

        mst.G <- mst.G + new.mst.G
    }
    mst.G[which(mst.G > 0)] = 1
    mst.G <- matrix2dgCMatrix(mst.G)
    return(mst.G)
}
