\name{HGC.PlotParameter}
\alias{HGC.PlotParameter}
\title{
Visualizing the Parameter Records during Clustering
}
\description{
The function visualizes the parameter output from \code{HGC.parameter}.
}
\usage{
HGC.PlotParameter(record, parameter)
}
\arguments{
    \item{record}{
    the input record matrix of parameters from \code{HGC.parameter}.
    }
    \item{parameter}{
    a string with alternatives "CL" or "ANN". Choose "CL" to plot the
    chain lengths and "ANN" to plot the average neighbor number.
    }
}
\details{
The chain length(CL) and average neighbor number(ANN) are key factors
related
with the time complexity of clustering by \code{HGC}. The function
provides the
visualization of them.
}
\value{
The function will return 1 if the dendrogram is sucessfully drawn.
}
\examples{
data(Pollen)
Pollen.PCs <- Pollen[["PCs"]]
Pollen.SNN <- SNN.Construction(Pollen.PCs)
Pollen.ParameterRecord <- HGC.parameter(G = Pollen.SNN)
HGC.PlotParameter(Pollen.ParameterRecord, parameter = "CL")
HGC.PlotParameter(Pollen.ParameterRecord, parameter = "ANN")
}
