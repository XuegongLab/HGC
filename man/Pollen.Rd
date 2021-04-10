\docType{data}
\name{Pollen}
\alias{Pollen}
\title{Embeddings of the Pollen datasets in the principal component space.}
\format{
An object of class \code{list}.

The \code{list} contains three elements: First, a matrix with
301 rows and 25 columns, saved as \code{Pollen[["PCs"]]}. The
rows are cells and columns are principal components. Second,
the label in tissue level, saved as a vector in
\code{Pollen[["Tissue"]]}.
Third, the label in cell line level, saved as a vector in
\code{Pollen[["CellLine"]]}.
}
\source{
\url{https://www.nature.com/articles/nbt.2967}
}
\description{
The dataset is the low dimensional principal components and
labels of cells from the Pollen dataset. The 301 cells are from
11 different cell lines and are classified into 4 tissues.
}
\keyword{datasets}
