# Rcpp added to ...
# ExGaussian stop-signal model with 1 stop & 2 go accumulators
rm(list=ls())
source ("dmc/dmc.R")
load_model ("EXG-SS","exgSS.R")


model <- model.dmc(
    # SS stands for trial type (GO or Stop-signal [SS]):
    factors=list(S=c("s1","s2"),SS=c("GO","SS")),
    # NR stands for "No response", i.e., go omission & successful inhibitions:
    responses=c("NR","r1","r2"),
    # Match scores correct responses for each GO stimulus as usual, and also scores
    # the "correct" stimulus corresponding to an NR response, but the latter has
    # no effect (it just avoids a standard check making sure that each response is scored):
    match.map=list(M=list(s1="r1",s2="r2",s1="NR")),
    p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
    type="exgss")

  # This parameter setting will produce about 25% (relatively slow) errors:
  p.vector  <- c(mu.true=.5,mu.false=.60,muS=.2,
                 sigma.true=.05,sigma.false=.030,sigmaS=.03,
                 tau.true=.08,tau.false=.04,tauS=.05,
                 tf=qnorm(.1),gf=qnorm(.1))



n <- c(375,375,125,125)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25)),model)

# # Plot go RT distribution;
# # P(NA) = Probability of go omission;
# # Accuracy is computed as correct go/(all go - go omissions):
# correct <- as.numeric(data$S)==(as.numeric(data$R)-1)
# layout(1)
# plot.cell.density(data[data$SS=="GO",],
#                   C=correct[data$SS=="GO"],
#                   xlim=c(0,5),ymax=5,main="Go RTs")
#
# # Overall accuracy (i.e., all trials-(errors + omission)
# tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]
# #  GO
# # 0.676
#
# # Show the different SSDs:
# sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS)
# # "0.15" "0.2"  "0.25" "0.3"  "0.35" "0.4"  "0.45"
#
# # Show the number of trials for each SSD:
# Ns = tapply(data$RT,data$SSD,length)
# Ns
# # 0.15  0.2 0.25  0.3 0.35  0.4 0.45  Inf
# #    1   10   43   79   77   35    5  750
#
# # Show response rate:
# tapply(!is.na(data$RT),data[,c("SS")],mean)
# #   GO    SS
# # 0.896 0.500
#
# # Response rate broken down by SSD & corresponding inhibition function:
# tapply(!is.na(data$RT),data[,c("SS","SSD")],mean)
# layout(1)
# plot_SS_if.dmc(data)  # P(Respond) increases as a function of SSD, as it should
# Ns
#
# # Median signal-respond RT per SSD:
# tapply(data$RT,data$SSD,median,na.rm=TRUE)
# layout(1)
# plot_SS_srrt.dmc(data) # Median SRRT increases as a function of SSD, as it should
# Nr = tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
# Nr
#
# # Signal-respond RTs should be faster than go RTs:
# hist(data$RT[data$SS=="GO"],breaks="fd",main="SRRT vs. Go RT",freq=F,ylim=c(0,8))
# lines(density(data$RT[data$SS=="SS"],na.rm=T),col="red",lwd=2)

### Testing, didnt try fitting here ----


source("dmc/models/EXG-SS/dists.R")
rlike  <- likelihood.dmc(p.vector,data)
Rcpp::sourceCpp(file = "dmc/models/EXG-SS/dists.cpp")
clike  <- likelihood.dmc(p.vector,data)
summary(rlike-clike)


# Truncated normal priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1[10:11] <- 3
p.prior <- prior.p.dmc(
  dists = rep("tnorm",length(p1)),p1=p.vector,p2=p1,
  lower=c(rep(0,length(p1)-2),-12,-12),upper=c(rep(2,3),rep(.5,6),rep(12,2))
)
par(mfcol=c(2,6)); for (i in names(p.prior)) plot.prior(i,p.prior)

# Start sampling
samples <- samples.dmc(nmc=100,p.prior,data)

source("dmc/models/EXG-SS/dists.R")
system.time({
  samplesR3 <- run.dmc(samples,report=1,p.migrate=.05)
})
#    user  system elapsed
# 222.620   9.245 276.259

Rcpp::sourceCpp(file = "dmc/models/EXG-SS/dists.cpp")
system.time({
  samplesC3 <- run.dmc(samples,report=1,p.migrate=.05)
})
 #   user  system elapsed
 # 90.062   2.835 100.005

