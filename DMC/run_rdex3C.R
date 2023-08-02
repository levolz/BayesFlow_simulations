rm(list=ls())
source ("dmc/dmc.R")
load_model ("WALD-SSEXG","waldSSexg.R")
Rcpp::sourceCpp(file = "dmc/models/WALD-SSEXG/dists.cpp")
load("rdex3C.RData")
system.time({
  samplesC3c <- RUN.dmc(samplesC3)
})
save(data,samples,samplesC3,samplesC3c,file="rdex3C.RData")
