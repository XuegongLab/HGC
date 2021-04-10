gg_color_ggplot <- function(n){
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

HGC.PlotDendrogram <- function(tree, k = 5, plot.label = FALSE, labels){
    checkDataStructure(tree, "tree", "hclust")
    ncell <- length(tree$order)
    if(ncell == 0){
        stop("There is no element in the tree.")
    }
    checkInteger(k, "k")
    if(k > ncell){
        stop("k should be smaller than cell number.")
    }
    if(plot.label){
        if(nrow(labels) != ncell){
            stop("The labels do not match the data.")
        }
    }
    clus.dendrogram <- stats::as.dendrogram(tree)
    dend <- clus.dendrogram %>%
            dendextend::set("branches_k_color", k=k) %>%
            dendextend::set("branches_lwd", 1.2) %>%
            dendextend::set("labels_colors", "white")
    if(plot.label == FALSE){
        plot(dend)
        return(1)
    }else{
        n_label <- ncol(labels)
        n_diff_labels <- rep(0, n_label)
        for(i in seq_len(n_label)){
            temp_n_label <- length(unique(labels[,i]))
            n_diff_labels[i] <- temp_n_label
        }
        T.colourbar <- matrix(nrow = ncell, ncol = n_label)
        colnames(T.colourbar) = colnames(labels)
        for(i in seq_len(n_label)){
            temp_label <- labels[,i]
            temp_u_label <- unique(labels[,i])
            temp_n_label <- n_diff_labels[i]
            temp_colourbar <- gg_color_ggplot(temp_n_label)
            for(j in seq_len(temp_n_label)){
                T.colourbar[which(temp_label ==
                                temp_u_label[j]),i] = temp_colourbar[j]
            }
        }
        plot(dend)
        dendextend::colored_bars(colors = T.colourbar,dend = dend)
        return(1)
    }
}

HGC.PlotARIs <- function(tree, k.min = 2, k.max = 15,
                        labels, return.ARI = FALSE){
    n_label <- ncol(labels)
    checkDataStructure(tree, "tree", "hclust")
    ncell <- length(tree$order)
    if(ncell == 0){
        stop("There is no element in the tree.")
    }
    if(k.min < 1 || round(k.min) != k.min || k.max < 1 ||
        round(k.max) != k.max){
        stop("k.min and k.max should be positive integer.")
    }
    if(k.min > ncell || k.max > ncell){
        stop("k.min and k.max should be smaller than cell number.")
    }
    if(k.min > k.max){
        stop("k.min should not be larger than k.max.")
    }
    if(nrow(labels) != ncell){
        stop("The labels do not match the data.")
    }
    hc.ari.mat1 <- matrix(nrow = k.max - k.min + 1, ncol = n_label)
    colnames(hc.ari.mat1) <- colnames(labels)
    for(j in seq_len(n_label)){
        temp_label <- labels[,j]
        for(i in seq(k.min, k.max)){
            hc.ari.mat1[i-k.min+1,j] = mclust::adjustedRandIndex(
                                        temp_label, cutree(tree, i))
        }
    }
    # HC line chart
    plot.df2 <- data.frame(k = factor(rep(seq(k.min, k.max), times = n_label)),
                            ARI = as.vector(hc.ari.mat1),
                            label = c(rep(colnames(hc.ari.mat1), each =
                                        k.max - k.min + 1)))

    p1 <- ggplot2::ggplot(plot.df2, ggplot2::aes_string(x="k", y="ARI",
                                        group="label", colour="label")) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size=1) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    plot(p1)
    if(return.ARI){
        return(hc.ari.mat1)
    }else{
        return(1)
    }
}

HGC.PlotParameter <- function(record, parameter = "CL"){
    n_iter <- ncol(record)
    checkDataStructure(record, "record", "matrix")
    if(nrow(record) != 2){
        stop("Wrong input data structure for record.")
    }
    if(rownames(record)[1] != "ChainLength" ||
        rownames(record)[2] != "AvgNeighNum"){
        stop("Wrong input data structure for record.")
    }
    if(parameter == "CL"){
        plot.df <- data.frame(Iter = seq_len(n_iter), ChainLength = record[1,])

        p1 <- ggplot2::ggplot(plot.df,
                            ggplot2::aes_string(x="Iter", y="ChainLength")) +
            ggplot2::geom_point(colour = "#00BFC4", size = 1.0)
        p2_hist <- ggplot2::ggplot(plot.df,
                            ggplot2::aes_string(x="ChainLength")) +
            ggplot2::geom_histogram(binwidth = 1,
                                    colour = "black",
                                    fill = "#00BF7D") +
            ggplot2::coord_flip() +
            ggplot2::theme_void()
        plot(patchwork::wrap_plots(p1, p2_hist,
                                    nrow = 1, widths = c(1, 0.5)))
        return(1)
    }else if(parameter == "ANN"){
        plot.df <- data.frame(Iter = seq_len(n_iter), ANN = record[2,])

        p1 <- ggplot2::ggplot(plot.df,
                            ggplot2::aes_string(x="Iter", y="ANN")) +
            ggplot2::geom_point(colour = "#F8766D", size = 1.0)
        p2_hist <- ggplot2::ggplot(plot.df,
                            ggplot2::aes_string(x="ANN")) +
            ggplot2::geom_histogram(binwidth = 1,
                                    colour = "black",
                                    fill = "#A3A500") +
            ggplot2::coord_flip() +
            ggplot2::theme_void()
        plot(patchwork::wrap_plots(p1, p2_hist,
                                    nrow = 1, widths = c(1, 0.5)))
        return(1)
    }else{
        stop("Wrong parameter name.")
    }
}
