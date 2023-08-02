rm(list=ls())
source ("dmc/dmc.R")
load_model ("WALD-SSEXG","waldSSexg.R")
load("rdex3R.RData")
system.time({
  samplesR3c <- RUN.dmc(samplesR3)
})
save(data,samples,samplesR3,samplesR3c,file="rdex3R.RData")
