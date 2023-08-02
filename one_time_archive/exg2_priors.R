rm(list = ls())

print(load("corey_data/corey_noError_hier.RData"))
# # corey is individual fits, h..1-4 is hierarchical, 4 is converged
# [1] "corey"   "hcorey"  "hcorey2" "hcorey3" "hcorey4"

thetas = lapply(corey, \(x) {
  matrix(aperm(x$theta, c(1,3,2)), ncol=8,
               dimnames = list(NULL, dimnames(x$theta)[[2]]))
})

theta = do.call(rbind, thetas)

corey_fits = list()
for (i in dimnames(theta)[[2]]) {
  print(i)
  # natural
  corey_fits[["natural"]][["full"]][[i]] = range(theta[, i])
  corey_fits[["natural"]][["95"]][[i]] = quantile(theta[, i], 
    probs = c(.025, .975), na.rm = TRUE)
  corey_fits[["natural"]][["mean"]][[i]] = mean(theta[, i], na.rm = TRUE)
  corey_fits[["natural"]][["sd"]][[i]] = sd(theta[, i], na.rm = TRUE)
  # log
  corey_fits[["log"]][["full"]][[i]] = range(log(theta[, i]))
  corey_fits[["log"]][["95"]][[i]] = quantile(log(theta[, i]), 
    probs = c(.025, .975), na.rm = TRUE)
  corey_fits[["log"]][["mean"]][[i]] = mean(log(theta[, i]), na.rm = TRUE)
  corey_fits[["log"]][["sd"]][[i]] = sd(log(theta[, i]), na.rm = TRUE)
}

# narural scale summary
corey_fits[["natural"]][["means"]] = c(
  corey_fits$natural$mean$mu.true, corey_fits$natural$mean$sigma.true,
  corey_fits$natural$mean$tau.true, corey_fits$natural$mean$muS,
  corey_fits$natural$mean$sigmaS, corey_fits$natural$mean$tauS
)
corey_fits[["natural"]][["sds"]] = c(
  corey_fits$natural$sd$mu.true, corey_fits$natural$sd$sigma.true,
  corey_fits$natural$sd$tau.true, corey_fits$natural$sd$muS,
  corey_fits$natural$sd$sigmaS, corey_fits$natural$sd$tauS
)
corey_fits[["natural"]][["cov"]] = cov(theta)[1:6, 1:6]  # mvn cov matrix


# log scale summary
corey_fits[["log"]][["means"]] = c(
  corey_fits$log$mean$mu.true, corey_fits$log$mean$sigma.true,
  corey_fits$log$mean$tau.true, corey_fits$log$mean$muS,
  corey_fits$log$mean$sigmaS, corey_fits$log$mean$tauS
)
corey_fits[["log"]][["sds"]] = c(
  corey_fits$log$sd$mu.true, corey_fits$log$sd$sigma.true,
  corey_fits$log$sd$tau.true, corey_fits$log$sd$muS, 
  corey_fits$log$sd$sigmaS, corey_fits$log$sd$tauS
)
corey_fits[["log"]][["cov"]] = cov(log(theta))[1:6, 1:6]  # mvn cov matrix


# save to RData for import to Python
save(corey_fits, file = "corey_data/corey_indiv-fits_extracted.RData")