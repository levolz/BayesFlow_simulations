rm(list=ls())
source("dmc/dmc.R")
load_model("EXG-SS", "exgSS.R")
Rcpp::sourceCpp(file = "dmc/models/EXG-SS/dists.cpp")
load("exg2C.RData")
system.time({
  samplesC2c <- RUN.dmc(samplesC2)
})
save(data, samples, samplesC2, samplesC2c, file = "exg2C.RData")