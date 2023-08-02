rm(list=ls())
source ("dmc/dmc.R")

load_model ("WALD-SSEXG","waldSSt50exgABCD_probit.R") 

model <- model.dmc(type="waldss",              
    factors=list(S=c("blue","orange"),SS=c("GO","SS")),
    responses = c("NR","BLUE","ORANGE"),
    match.map=list(M=list(blue="BLUE",orange="ORANGE",orange="NR")),
    p.map=list(A="1",B="1",v=c("M"), t0="1", gf="1",
               mu="1",sigma="1",tau="1",tf="1",ts="1",
               vT="1",vF="1",v0="1",aT="1",aF="1"), # Initial common rate at SSD=0, exponential 
                                      # rates of change with SSD for TRUE and FALSE
    constants = c(A=0,ts=0,v.true=Inf,v.false=-Inf)) 
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "B"     "t0"    "gf"    "mu"    "sigma" "tau"   "tf"   
#  [8] "vT"    "vF"    "v0"    "aT"    "aF"   
# 
# Constants are (see attr(,"constants") ):
#       A      ts  v.true v.false 
#       0       0     Inf    -Inf 
# 
# Model type = waldss  

p.vector  <- c(vT=2.5,vF=.5,B=1,mu=0.3,sigma=0.05,tau=0.1,
               t0=.2,gf=qnorm(.1),tf=qnorm(.1),v0=1,aT=1,aF=1)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)
dat <- simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25))
data <- data.model.dmc(dat,model) 
# likelihood.dmc(p.vector,data)

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,6))
profile.dmc("vT",2,3,p.vector,data)
profile.dmc("vF",.25,0.75,p.vector,data)
profile.dmc("mu",.2,.4,p.vector,data)
profile.dmc("sigma",.01,.075,p.vector,data)
profile.dmc("tau",.05,.15,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("aT",.5,1.5,p.vector,data)
profile.dmc("aF",.5,1.5,p.vector,data)
profile.dmc("v0",.5,1.5,p.vector,data)
profile.dmc("tf",qnorm(.05),qnorm(.15),p.vector,data)
profile.dmc("gf",qnorm(.05),qnorm(.15),p.vector,data)


