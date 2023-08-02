rm(list = ls())

corey_folder = "corey_data"


print(load(paste(corey_folder, "corey_noError_hier.RData", sep = "/")))
# # corey is individual fits, h..1-4 is hierarchical, 4 is converged
# [1] "corey"   "hcorey"  "hcorey2" "hcorey3" "hcorey4"

# Pull out the data, 25% stop (128 trials, 96 go, 32 stop)
data <- lapply(corey,function(x)x$data)
ns <- unlist(lapply(data,function(x)dim(x)[1]))
dat <- cbind(subjects=rep(names(corey),ns),do.call(rbind,data))
row.names(dat) <- NULL
# Make in to EMC format
names(dat)[5] <- "rt"
dat <- dat[,c(1,2,4,6,5)]
levels(dat$R) <- c("stop","left","right")
levels(dat$S) <- c("left","right")

# # Some checks
# # Only correct responses kept
# tmp <- dat[!is.na(dat$rt),]; all(as.character(tmp$S)==as.character(tmp$R))
# # Ignoring response looks fairly justified.
# # right ~ 7ms faster on average, but fairly identical in distribution
# tapply(tmp$rt,tmp$S,mean)
# qqplot(tmp$rt[tmp$S=="left"],tmp$rt[tmp$S=="right"]);abline(a=0,b=1,col="red",lwd=2)


# Pull out list of estimates
thetas = lapply(corey,function(x) {
  matrix(aperm(x$theta, c(1,3,2)), ncol=8,
               dimnames = list(NULL, dimnames(x$theta)[[2]]))
})

# Put into a data frame
theta = do.call(rbind, thetas)

# Box-Cox ladder of powers
boxcox <- function(te,p) if (te==0) log(p) else p^te

untrans <- function(p,te) {
	do_log <- te==0
	out <- p
	out[!do_log] <- p[!do_log]^(1/te[!do_log])
	out[do_log] <- exp(p[do_log])
	do_probit <- names(out) %in% c("tf","gf")
	out[do_probit] <- pnorm(out[do_probit])
	out
}

# maximize W normality statistic (max=1)
obj <- function(te,p) {
	shapiro.test(quantile(boxcox(te,p),probs=c(1:999)/1000))$statistic
}

get_trans <- function(theta) {
  trans <- setNames(vector(mode="list",length=dim(theta)[2]),dimnames(theta)[[2]]) 
  for (i in dimnames(theta)[[2]]) { 
	  if (!(i %in% c("tf","gf"))) {
	    trans[[i]] <- optimize(obj,c(-4,4),p=theta[,i],maximum=TRUE)
  	} else trans[[i]] <- list(maximum=1,
	    objective=shapiro.test(quantile(theta[,i],probs=c(1:999)/1000))$statistic)
  }
  trans
}


# Plot transforms
plot_trans <- function(trans,theta) {
	theta_trans <- theta
  par(mfrow=c(2,4))
  for (i in dimnames(theta)[[2]]) { 
	  theta_trans[,i] <- boxcox(trans[[i]]$maximum,theta[,i])
  	x <- hist(theta_trans[,i],breaks="fd",freq=FALSE,
	    main=paste(i,"BC =",round(trans[[i]]$maximum,2),"W =",round(trans[[i]]$objective,3)))
  	lines(x$mids,dnorm(x$mids,mean(theta_trans[,i]),sd(theta_trans[,i])),col="red",lwd=2)
  }
  attr(theta_trans,"trans") <- trans
  invisible(theta_trans)
}

update_trans <- function(te,trans,theta) {
	for (i in names(te)) {
	  trans[[i]] <- list(maximum=te[i],
	      objective=shapiro.test(quantile(boxcox(te[i],theta[,i]),probs=c(1:999)/1000))$statistic)
	}
	trans
}

lo_windsor <- function(theta,lo) {
	for (i in names(lo)) {
		bad <- theta[,i]< lo[i]
	  theta[bad,i] <- lo[i]	
	  bad <- 100*mean(bad)
	  names(bad) <- i
	  print(bad)
	}
	theta
}

hi_windsor <- function(theta,hi) {
	for (i in names(hi)) {
		bad <- theta[,i]> hi[i]
	  theta[bad,i] <- hi[i]	
	  bad <- 100*mean(bad)
	  names(bad) <- i
	  print(bad)
	}
	theta
}


# Get optimal transform	
trans <- get_trans(theta)
# Get no transform (te=1)
trans1 <- update_trans(setNames(rep(1,dim(theta)[2]),dimnames(theta)[[2]]),trans,theta)

# Look at results
theta0 <- theta
theta0[,c("tf","gf")] <- pnorm(theta0[,c("tf","gf")])
round(apply(theta0,2,range),3)
round(apply(theta0,2,quantile,probs=c(.01,.99)),3)


# Plot untransformed and optimal
plot_trans(trans1,theta)
plot_trans(trans,theta)
theta_trans <- plot_trans(trans,theta)
# extreme values in mu.true, due to implausibly fast values, Windsorise at .2s
windsor <- c(mu.true=.2)
theta_windsor <- lo_windsor(theta,windsor)
trans_windsor <- get_trans(theta_windsor)
# Little effect on optimal transform
theta_trans_windsor <- plot_trans(trans_windsor,theta_windsor)

te_final <- setNames(c(-1,0,.5,.5,0,0,1,1),dimnames(theta)[[2]])
trans_final <- update_trans(te_final,trans,theta)

theta_trans_final <- plot_trans(trans_final,theta)
theta_trans_windsor_final <- plot_trans(trans_final,theta_windsor)


library(MASS)
library(tmvnsim)

untrans_mat <- function(p,te) {
	do_log <- te==0
	out <- p
	do_probit <- dimnames(out)[[2]] %in% c("tf","gf")
	do_power <- !do_log & !do_probit
	
	out[,do_log] <- exp(p[,do_log])
	# Leave as probit
	# out[,do_probit] <- pnorm(out[,do_probit])
	for (i in dimnames(out)[[2]][do_power])
		out[,i] <- p[,i]^(1/te[i])
	out
}

theta_final <- theta_trans_windsor_final
mu <- apply(theta_final,2,mean)
Sigma <- cov(theta_final)
round(cov2cor(Sigma),2)

# This produces rare unreasonably large samples for mu/sigma/tau
samps <- mvrnorm(1e6,mu,Sigma)
# Truncate samles above at 1s (the timeout of the experiment)
# Also bound tf and gf above by 0.5 as otherwise little useful data
tmv_lower = c(1,-Inf, rep(0,2), rep(-Inf,4))
tmv_upper = c(Inf,0,rep(1,2),rep(0,2),qnorm(.5),qnorm(.5))
samps <- tmvnsim(1e6,length(mu),means=mu,sigma=Sigma,
  lower=tmv_lower,upper=tmv_upper)$samp
dimnames(samps)[[2]] <- names(mu)


theta_sim <- untrans_mat(samps,te_final)
trans_sim <- update_trans(te_final,trans,theta_sim)
# Look at the raw samples
plot_trans(trans1,theta_sim)
# Confirm the transformed samples are marginally (sometiems turncated) gaussian
plot_trans(trans_sim,theta_sim)

# Look at real ranges
tmp <- theta_sim; tmp[,c("tf","gf")] <- pnorm(tmp[,c("tf","gf")])
round(apply(tmp,2,range),3)

# Final self-contained simulator

corey_pars <- function(n,mu,Sigma,te=NULL,
	lower=c(1,rep(-Inf,7)),upper=c(Inf,0,rep(1,2),rep(0,2),qnorm(.5),qnorm(.5))) {

	muntrans <- function(p,te) {
  	do_log <- te==0
	  out <- p
  	do_probit <- dimnames(out)[[2]] %in% c("tf","gf")
	  do_power <- !do_log & !do_probit
  	out[,do_log] <- exp(p[,do_log])
  	out[,do_probit] <- pnorm(out[,do_probit])
	  for (i in dimnames(out)[[2]][do_power])
		  out[,i] <- p[,i]^(1/te[i])
  	out
}

	samps <- tmvnsim(n,length(mu),means=mu,sigma=Sigma,lower=lower,upper=upper)$samp
  dimnames(samps)[[2]] <- names(mu)
  if (!is.null(te)) muntrans(samps,te) else samps
}

# Make transformed
cp <- corey_pars(1e6,mu,Sigma)
#Check
plot_trans(trans1,cp)
# Make natural
cpt <- corey_pars(1e6,mu,Sigma,te=te_final)
# Check
plot_trans(trans1,cpt)

save(mu,Sigma,te_final,trans1,corey_pars,plot_trans,boxcox,tmv_lower,tmv_upper,
	 file=paste(corey_folder, "corey_simulator.RData", sep = "/"))