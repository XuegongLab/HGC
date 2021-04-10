## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE------------------
library(knitr)
opts_chunk$set(fig.align = 'center', 
                fig.width = 4.5, fig.height = 3, dev = 'png')

## ----Bioconductor install, eval = FALSE---------------------------------------
#  if (!requireNamespace("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install("HGC")

## ----Github install, eval = FALSE---------------------------------------------
#  if(!require(devtools))
#      install.packages("devtools")
#  devtools::install_github("stevenhuakui/HGC")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(HGC)

data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
Pollen.Label.Tissue <- Pollen[["Tissue"]]
Pollen.Label.CellLine <- Pollen[["CellLine"]]

dim(Pollen.PCs)
table(Pollen.Label.Tissue)
table(Pollen.Label.CellLine)

## -----------------------------------------------------------------------------
Pollen.SNN <- SNN.Construction(mat = Pollen.PCs, k = 25, threshold = 0.15)
Pollen.ClusteringTree <- HGC.dendrogram(G = Pollen.SNN)

## -----------------------------------------------------------------------------
cluster.k5 <- cutree(Pollen.ClusteringTree, k = 5)

## ---- fig.height = 4.5--------------------------------------------------------
Pollen.ClusteringTree$height = log(Pollen.ClusteringTree$height + 1)
Pollen.ClusteringTree$height = log(Pollen.ClusteringTree$height + 1)

HGC.PlotDendrogram(tree = Pollen.ClusteringTree,
                    k = 5, plot.label = FALSE)

## ---- fig.height = 4.5--------------------------------------------------------
Pollen.labels <- data.frame(Tissue = Pollen.Label.Tissue,
                            CellLine = Pollen.Label.CellLine)
HGC.PlotDendrogram(tree = Pollen.ClusteringTree,
                    k = 5, plot.label = TRUE, 
                    labels = Pollen.labels)

## -----------------------------------------------------------------------------
ARI.mat <- HGC.PlotARIs(tree = Pollen.ClusteringTree,
                        labels = Pollen.labels)

## -----------------------------------------------------------------------------
Pollen.ParameterRecord <- HGC.parameter(G = Pollen.SNN)

HGC.PlotParameter(Pollen.ParameterRecord, parameter = "CL")
HGC.PlotParameter(Pollen.ParameterRecord, parameter = "ANN")

