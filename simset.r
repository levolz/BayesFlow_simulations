rm(list = ls())

source("dmc/dmc/dmc.R")
load_model("exg-ss", "exgSS.R", locale = "dmc/dmc/models")
Rcpp::sourceCpp(file = "dmc/dmc/models/EXG-SS/dists.cpp")

#### constants ####
PARAM_NAMES <- c("mu", "sigma", "tau", "muS", "sigmaS", "tauS")

CONSTANTS <- list(
  mu.false = 1e-6, # box-cox transformed
  sigma.false = -9, # box-cox transformed
  tau.false = 1e-2, # box-cox transformed
  tf = 0, gf = 0, # we simulate with no failures
  lower = c(1, -1e6, 1e-6, 1e-6, -1e6, -1e6),
  upper = c(1e6, 0, 1, 1, 1e-6, 0),
  transform_vec = c(-1, 0, .5, .5, 0, 0)
)


#### Helpers ####
untrans <- function(p, te) {
  do_log <- (te == 0)
  out <- p
  out[, !do_log] <- t(t(p[, !do_log])^(1 / te[!do_log]))
  out[, do_log] <- exp(p[, do_log])
  do_probit <- names(out) %in% c("tf", "gf")
  out[, do_probit] <- pnorm(out[, do_probit])
  out
}

duntrans <- function(p, te = transform_vec) {
  do_log <- names(te)[te == 0]
  do_pow <- names(te)[te != 0]
  out <- p
  for (i in do_log) {
    out[, i] <- exp(p[, i])
  }
  for (i in do_pow) {
    out[, i] <- p[, i]^(1 / te[i])
  }
  out
}

untrans_fun <- function(te) {
  do_log <- (te == 0)
  out <- te
  out[!do_log] <- \(x) x^(1 / te[!do_log])
  out[do_log] <- exp
  do_probit <- names(out) %in% c("tf", "gf")
  out[do_probit] <- pnorm(out)
}

unpack_bf2dmc <- function(data, s = NULL) {
  # unpacks data from bf2dmc format to dmc format
  # data: data.frame with columns "rt", "ssd"
  # returns: data.frame with columns "RT", "SSD", "s", "S", "SS", "R"
  colnames(data) <- c("RT", "SSD")

  if (!is.null(s)) data$s <- rep(s, nrow(data))
  data$S <- c("s1", "s2")
  data$SS <- ifelse(data$SSD == -1.0, "GO", "SS") # CHANGE TO -1.0
  data$SSD[data$SSD == -1.0] <- Inf # CHANGE TO -1.0

  data$RT[data$RT == -1.0] <- NA # CHANGE TO -1.0
  data$R <- c("r1", "r2")
  data$R[is.na(data$RT)] <- "NR"

  data$S <- factor(data$S, levels = c("s1", "s2"))
  data$SS <- factor(data$SS, levels = c("GO", "SS"))
  data$R <- factor(data$R, levels = c("NR", "r1", "r2"))

  return(data)
}

unpack_dmc2bf <- function(data) {
  # unpacks data from dmc format to bf format
  # data: data.frame with columns "S", "SS", "R", "RT", "SSD"
  # returns: data.frame with columns "RT", "SSD"
  data$RT[is.na(data$RT)] <- -1.0
  data$SSD[data$SSD == Inf] <- -1.0
  return(data[, c("RT", "SSD")])
}

extract_thetas <- function(fits, ci = 0.95) {
  mean <- apply(fits$theta, MARGIN = 2, mean)
  sd <- apply(fits$theta, MARGIN = 2, sd)
  median <- apply(fits$theta, MARGIN = 2, median)
  low <- apply(fits$theta, MARGIN = 2, quantile, probs = .5 * (1 - ci))
  high <- apply(fits$theta, MARGIN = 2, quantile, probs = .5 * (1 + ci))
  return(data.frame(
    mean = mean, sd = sd,
    median = median, low = low, high = high
  ))
}

transform_theta_fits <- function(
    fit,
    trans = transform_vec,
    summary = "median") {
  thetas <- data.frame()
  for (i in seq_len(dim(fit$theta)[1])) {
    thetas <- rbind(thetas, untrans(t(fit$theta[i, , ]), trans))
  }
  out <- list()
  if ("median" %in% summary) {
    out$median <- apply(thetas, 2, median)
  }
  if ("mean" %in% summary) {
    out$mean <- colMeans(thetas)
  }
  if ("sd" %in% summary) {
    out$sd <- apply(thetas, 2, sd)
  }
  if ("mad" %in% summary) {
    out$mad <- apply(thetas, 2, mad)
  }
  return(out)
}

get_corey_data_to_bf <- function(data, n = length(data)) {
  dmc_data <- lapply(data, \(x) x$data)
  bf_data <- lapply(
    dmc_data, unpack_dmc2bf
  )
  # write each entry to csv
  for (i in names(bf_data)) {
    write.csv(bf_data[[i]], paste0("corey_data/bf_new/corey_", i, ".csv"), row.names = F)
  }
}
