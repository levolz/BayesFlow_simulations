# RDEX3 with rcpp

rm(list=ls())
source ("dmc/dmc.R")
load_model ("WALD-SSEXG","waldSSexg.R")

# With gf and tf, A=0
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters
             t0="1",                           # GO non-decision
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             minEXG="1",                       # min for stop distribution
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  # constants = c(A=0,ts=0)) # No start point noise or TRIAL effects
  constants = c(A=0,ts=0,minEXG=0)) # No start point noise or TRIAL effects

p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,tf=.1,gf=.1)


n <- c(375,375,125,125)
dat <- simulate.dmc(p.vector,model,staircase=.05,
                    n=n,SSD=c(Inf,Inf,.25,.25))
data <- data.model.dmc(dat,model)

# check likelihoods are the same ----

load_model ("WALD-SSEXG","waldSSexg.R")
rlike  <- likelihood.dmc(p.vector,data)
Rcpp::sourceCpp(file = "dmc/models/WALD-SSEXG/dists.cpp")
clike  <- likelihood.dmc(p.vector,data)
summary(rlike-clike)

# Truncated normal priors:
p1 <- p.vector
p2 <- rep(1,length(p1))
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p1,
  lower=c(rep(0,length(p1)-2),-12,-12),upper=c(rep(5,3),rep(1,4),rep(12,2))
)
par(mfcol=c(2,5)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)

load_model ("WALD-SSEXG","waldSSexg.R")
system.time({
  samplesR3 <- run.dmc(samples,report=1,p.migrate=.05)
})
 #   user  system elapsed
 # 54.426   0.491  54.989
save(data,samples,samplesR3,file="rdex3R.RData")

Rcpp::sourceCpp(file = "dmc/models/WALD-SSEXG/dists.cpp")
system.time({
  samplesC3 <- run.dmc(samples,report=1,p.migrate=.05)
})
 #   user  system elapsed
 # 20.922   0.171  21.086
save(data,samples,samplesC3,file="rdex3C.RData")


# # RUN.dmc in run_exg2.R
#    user  system elapsed
# 299.829   1.920 302.427
# #  and run_exg2C.R
#    user  system elapsed
# 116.535   1.136 118.213
# 2.56 speedup

print(load("rdex3R.RData"))
print(load("rdex3C.RData"))

par(mfrow=c(1,2))
plot.dmc(samplesR3c,pll.chain=TRUE)
plot.dmc(samplesC3c,pll.chain=TRUE)

# Parameter chains look like fat hairy caterpillars:
plot.dmc(samplesR3c,layout=c(2,5))
plot.dmc(samplesC3c,layout=c(2,5))

check.recovery.dmc(samplesR3c,p.vector)
#                v.true v.false    B    t0    mu sigma  tau    tf    gf
# True             2.50    0.50 1.00  0.20  0.30  0.05 0.10  0.10  0.10
# 2.5% Estimate    2.36    0.43 0.96  0.17  0.20  0.00 0.15 -0.73  0.08
# 50% Estimate     2.57    0.69 1.05  0.19  0.25  0.04 0.20 -0.58  0.09
# 97.5% Estimate   2.78    0.96 1.16  0.21  0.30  0.10 0.27 -0.43  0.11
# Median-True      0.07    0.19 0.05 -0.01 -0.05 -0.01 0.10 -0.68 -0.01
check.recovery.dmc(samplesC3c,p.vector)
#                v.true v.false    B    t0    mu sigma  tau    tf    gf
# True             2.50    0.50 1.00  0.20  0.30  0.05 0.10  0.10  0.10
# 2.5% Estimate    2.38    0.44 0.96  0.17  0.20  0.00 0.15 -0.73  0.08
# 50% Estimate     2.57    0.69 1.05  0.19  0.26  0.04 0.20 -0.58  0.09
# 97.5% Estimate   2.81    0.96 1.18  0.21  0.30  0.10 0.27 -0.44  0.11
# Median-True      0.07    0.19 0.05 -0.01 -0.04 -0.01 0.10 -0.68 -0.01




