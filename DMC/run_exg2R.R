rm(list = ls())
source("dmc/dmc.R")
load_model("EXG-SS", "exgSS.R")
load("exg2R.RData")
system.time({
  samplesR2c <- RUN.dmc(samplesR2)
})
save(data, samples, samplesR2, samplesR2c, file = "exg2R.RData")
