##################  DMC Stop Signal Model Lesson
#
rm(list=ls())
source ("dmc/dmc.R")
#
# This lesson illustrates the use of 
# 1) waldSSexg:  
# 2) waldSSesN: The same as (1) but with a variable number of accumulators
# 3) waldSSexg with slowing

# Use a pure Wald model (no start point variability, so A=0, AS=0).


# is.tf <- TRUE
# is.gf <- TRUE
# is.ts <- FALSE
# use.staircase <- TRUE

########################################################## 1) ----

load_model ("WALD-SSEXG","waldSSexg.R") 

# Basic Model
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,ts=0,tf=0,gf=0)) # No start point noise or TRIAL effects

p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2)


# With gf and tf
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,ts=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "v.true"  "v.false" "B"       "mu"      "sigma"   "tau"     "t0"     
# [8] "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
#  A ts 
#  0  0 
# 
# Model type = waldss

# 
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,tf=.1,gf=.1)

# Full model (didnt test ts = trial slope as slow!)
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(ts=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "v.true"  "v.false" "B"       "mu"      "sigma"   "tau"     "t0"     
# [8] "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
#  A ts 
#  0  0 
# 
# Model type = waldss

# 
p.vector  <- c(v.true=2.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,tf=.1,gf=.1,A=.5)


# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
n=6
data <- data.model.dmc(simulate.dmc(p.vector,model,
                       n=n,SSD=c(Inf,Inf,.25,.25)),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 

# mean(data$RT[data$SS=="GO"])
# mean(data$RT[data$SS!="GO"],na.rm=TRUE)

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("mu",.2,.4,p.vector,data)
profile.dmc("sigma",.01,.075,p.vector,data)
profile.dmc("tau",.05,.15,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)

profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)

profile.dmc("A",.25,.75,p.vector,data)

# recovery all good expect sigma consistently under-estiamted.

##################################  Mixed 2 and 3 accumulator example ----

load_model ("WALD-SSEXG","waldSSexgN.R") 

# Just test full model
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2","s3"),SS=c("GO","SS"), # Go stimulus, go and stop trials
               NC=c("a2","a3")),               # 2 choice vs. 3 choice case
  responses <- c("NR","r1","r2","r3"),         # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             N="NC",                           # Number of accumulators
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")), # NR mapping ignored
  constants = c(ts=0,N.a2=3,N.a3=4)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "v.true"  "v.false" "B"       "A"       "mu"      "sigma"   "tau"    
#  [8] "t0"      "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
#   ts N.a2 N.a3 
#    0    3    4 
# 
# Model type = waldss 


p.vector  <- c(v.true=3,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,tf=.1,gf=.1,A=.5)

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
# Full factorial 12 cells, but remove s3 in a2 for both GO and SS
SSD <- c(rep(Inf,3),rep(.25,3),rep(Inf,3),rep(.25,3)) # Per cell, two ignored
# SSD <- c(rep(Inf,2*6),rep(.25,2*6),rep(Inf,3*6),rep(.25,3*6)) # Per value

# Check likelihood
n <- c(6,6,0,6,6,0,6,6,6,6,6,6)  
data <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=SSD),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,0,2.5e3,2.5e3,0,7.5e3,7.5e3,7.5e3,2.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=SSD),model) 

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2.5,3.5,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("mu",.2,.4,p.vector,data)
profile.dmc("sigma",.025,.075,p.vector,data)
profile.dmc("tau",.05,.15,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)
profile.dmc("A",.25,.75,p.vector,data)


##################### 3) SLOWING EXAMPLE  ----

load_model ("WALD-SSEXG","waldSSexg.R") 

# Basic Model
model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),               # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             mu="1",sigma="1",tau="1",         # STOP accumulator parameters
             t0="1",                           # GO non-decision
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,tf=0,gf=0) # No start point noise or TRIAL effects
)

p.vector  <- c(v.true=1.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,ts=1) # B increases from 1 to 2 given max(TRIAL) = 1 

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
n=c(6,6,2,2)
TRIALS <- c(1:sum(n))
data <- data.model.dmc(simulate.dmc(p.vector,model,
                       n=n,SSD=c(Inf,Inf,.25,.25),TRIALS=TRIALS),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)/10
TRIALS <- c(1:sum(n))/sum(n) 
# This screws things up
TRIALS <- c(1:sum(n))[order(runif(sum(n)))]/sum(n)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25),TRIALS=TRIALS),model) 
# # data <- data[order(data$TRIAL),]
# tmp <- likelihood.dmc(p.vector,data)
# tail(tmp)
# tmp <- cbind(data,tmp); View(tmp)

# mean(data$RT[data$SS=="GO"])
# mean(data$RT[data$SS!="GO"],na.rm=TRUE)

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# c(v.true=1.5,v.false=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,t0=.2,ts=1)
# Check the profiles
par(mfrow=c(2,4))
profile.dmc("v.true",1,2,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("mu",.2,.4,p.vector,data)
profile.dmc("sigma",.01,.075,p.vector,data)
profile.dmc("tau",.05,.15,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("ts",.5,1.25,p.vector,data)

# profile.dmc("gf",.05,.15,p.vector,data)
# profile.dmc("tf",.05,.15,p.vector,data)
# profile.dmc("A",.25,.75,p.vector,data)
