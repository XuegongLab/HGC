set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

sim.SNN <- SNN.Construction(simdata)
sim.ClusteringTree <- HGC.dendrogram(G = sim.SNN)

# Good Input
sim.labels <- data.frame(Tissue = sample(c("A","B"), 100, replace = TRUE),
                            CellLine = sample(c("AA","BA", "BB"), 100, replace = TRUE))

test_that("The function could show the figure", {
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 5, plot.label = FALSE), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 10, plot.label = FALSE), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 50, plot.label = FALSE), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 2, plot.label = TRUE, labels = sim.labels), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 5, plot.label = TRUE, labels = sim.labels), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 10, plot.label = TRUE, labels = sim.labels), 1)
    expect_equal(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 50, plot.label = TRUE, labels = sim.labels), 1)
})

# Bad Input: tree
test_that("Wrong input for tree",{
    expect_error(HGC.PlotDendrogram(tree = 1, k = 2, plot.label = FALSE), "Wrong input data structure for tree.")
    expect_error(HGC.PlotDendrogram(tree = simdata, k = 2, plot.label = FALSE), "Wrong input data structure for tree.")
})

# Bad input: k
test_that("Wrong input for k",{
    expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = -1, plot.label = FALSE), "k should be positive integer.")
    expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 0, plot.label = FALSE), "k should be positive integer.")
    expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 1.01, plot.label = FALSE), "k should be positive integer.")
    expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 100000, plot.label = FALSE), "k should be smaller than cell number.")
})

# Bad input: labels
test_that("Wrong input for labels",{
  expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 2, plot.label = TRUE), "argument \"labels\" is missing, with no default")
  expect_error(HGC.PlotDendrogram(tree = sim.ClusteringTree, k = 2, plot.label = TRUE, labels = sim.labels[1:5,]), "The labels do not match the data.")
})
