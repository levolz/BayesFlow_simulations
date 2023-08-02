rm(list = ls())

source("simset.r")

library(parallel)
# library(parallelly)

load("corey_data/corey_simulator.RData")
transform_vec <- te_final[1:6]
names(transform_vec) <- PARAM_NAMES

file_name <- "uvn-dmc"
fitted_file <- paste("./validation_data/", file_name, "_fits.RData", sep = "")


#### load data ####
load("corey_data/corey_noError_hier.RData")
val_data <- lapply(corey, \(x) x$data)


#### fit dmc model ####

# EXG2 model data generation
model_exg2 <- model.dmc(
  type = "exgss",
  factors = list(S = c("s1", "s2"), SS = c("GO", "SS")),
  responses = c("NR", "r1", "r2"),
  p.map = list(
    mu = "M", sigma = "M", tau = "M",
    muS = "1", sigmaS = "1", tauS = "1",
    tf = "1", gf = "1"
  ),
  match.map = list(M = list(s1 = "r1", s2 = "r2", s1 = "NR")),
  constants = c(
    mu.false = CONSTANTS[["mu.false"]],
    sigma.false = CONSTANTS[["sigma.false"]],
    tau.false = CONSTANTS[["tau.false"]],
    tf = CONSTANTS[["tf"]],
    gf = CONSTANTS[["gf"]]
  )
)

p.vector <- mu[1:6]

# loads mu vector from fitted vals
# lower_bounds <- c(1, -1e6, 1e-6, 1e-6, -1e6, -1e6)
# upper_bounds <- c(1e6, 0, 1, 1, 1e-6, 0)

check.p.vector(p.vector, model_exg2)

prior_exg2 <- prior.p.dmc(
  dists = rep("tnorm", length(p.vector)),
  p1 = p.vector, p2 = sqrt(diag(Sigma)[1:6]),
  lower = CONSTANTS[["lower"]], upper = CONSTANTS[["upper"]]
)

if (TRUE) {
  data_model <- data.model.dmc(val_data[[1]], model_exg2)
  samples_recov <- samples.dmc(
    nmc = 120,
    data = data_model,
    p.prior = prior_exg2
  )
  system.time({
    fin <- RUN.dmc(samples_recov)
  })
  thetas <- data.frame()
  for (i in 1:dim(fin$theta)[1]) {
    thetas <- rbind(thetas, untrans(t(fin$theta[i, , ]), transform_vec))
  }
  colMeans(thetas)
}


# Run DMC
if (file.exists(fitted_file)) {
  load(fitted_file)
} else {
  # setup samples object
  data_models <- lapply(val_data, function(x) data.model.dmc(x, model_exg2))
  #
  samples_recov <- lapply(data_models, function(x) {
    samples.dmc(
      nmc = 120,
      data = x, p.prior = prior_exg2
    )
  })
  recov_fits <- lapply(samples_recov,
    FUN = function(x) {
      cat("run.dmc\n")
      RUN.dmc(x)
    }
  )
  save(recov_fits, file = fitted_file)
}


theta_summary <- lapply(recov_fits, transform_theta_fits, summary = c("median", "mad"))
theta_true <- lapply(corey, extract_thetas)

fitted_thetas <- t(sapply(theta_summary, \(x) x$median))
orig_fits <- t(sapply(theta_true, \(x) x$median[1:6]))

fitted_mad <- t(sapply(theta_summary, \(x) x$mad))

windows()
par(mfrow = c(2, 3))
for (i in seq_len(length(PARAM_NAMES))) {
  plot(orig_fits[, i], fitted_thetas[, i],
    xlab = "True", ylab = "Fitted",
    main = PARAM_NAMES[i]
  )
}

write.csv(
  fitted_thetas,
  paste("./corey_data/", file_name, "_params.csv",
    sep = ""
  ),
  row.names = FALSE, col.names = FALSE
)
write.csv(
  orig_fits,
  paste("./corey_data/", file_name, "_orig_params.csv",
    sep = ""
  ),
  row.names = FALSE, col.names = FALSE
)

write.csv(
  fitted_mad,
  paste("./corey_data/", file_name, "_mad.csv",
    sep = ""
  ),
  row.names = FALSE, col.names = FALSE
)
