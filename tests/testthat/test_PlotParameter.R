set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

sim.SNN <- SNN.Construction(simdata)
sim.ParameterRecord <- HGC.parameter(G = sim.SNN)

# Good Input
test_that("The function could show the figure", {
    expect_equal(HGC.PlotParameter(record = sim.ParameterRecord, parameter = "CL"), 1)
    expect_equal(HGC.PlotParameter(record = sim.ParameterRecord, parameter = "ANN"), 1)
})

# Bad Input: record
test_that("Wrong input for record",{
    expect_error(HGC.PlotParameter(record = 123, parameter = "CL"), "Wrong input data structure for record.")
    expect_error(HGC.PlotParameter(record = simdata, parameter = "CL"), "Wrong input data structure for record.")
    expect_error(HGC.PlotParameter(record = sim.SNN, parameter = "CL"), "Wrong input data structure for record.")
})

# Bad input: parameter
test_that("Wrong input for parameter",{
    expect_error(HGC.PlotParameter(record = sim.ParameterRecord, parameter = "XXX"), "Wrong parameter name.")
})
