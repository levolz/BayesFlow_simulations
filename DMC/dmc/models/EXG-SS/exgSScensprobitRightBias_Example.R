rm(list=ls())
# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")

# your model 
load_model ("EXG-SS","exgSScensprobit.R")
model <- model.dmc(
   factors=list(S=c("left","right"),SS=c("go","stop")),
   responses=c("NR","LEFT","RIGHT"),
   match.map=list(M=list(left="LEFT",right="RIGHT",left="NR")),
   p.map=list(mu="M",sigma="M",tau="M",muS="1",sigmaS="1",tauS="1",tf="1",UC="1"),
   constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,UC=1),
   type="exgss")

p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,
                 tf=qnorm(.1))

check.p.vector(p.vector,model)
n <- c(375,375,125,125)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25)),model)
likelihood.dmc(p.vector,data)

# muB,m sigBR and tauBR add something to the right accumulator

load_model ("EXG-SS","exgSScensprobitRbias.R")
model <- model.dmc(
   factors=list(S=c("left","right"),SS=c("go","stop")),
   responses=c("NR","LEFT","RIGHT"),
   match.map=list(M=list(left="LEFT",right="RIGHT",left="NR")),
   p.map=list(mu="M",sigma="M",tau="M",muS="1",sigmaS="1",tauS="1",tf="1",UC="1",
     muBR="1",sigBR="1",tauBR="1"),
   constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,UC=1),
   type="exgss")


p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,
                 tf=qnorm(.1),muBR=.1,sigBR=.1,tauBR=.1)


check.p.vector(p.vector,model)
n <- c(375,375,125,125)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25)),model)

likelihood.dmc(p.vector,data)

# To not do some of the three biases simply omit

load_model ("EXG-SS","exgSScensprobitRbias.R")
model <- model.dmc(
   factors=list(S=c("left","right"),SS=c("go","stop")),
   responses=c("NR","LEFT","RIGHT"),
   match.map=list(M=list(left="LEFT",right="RIGHT",left="NR")),
   p.map=list(mu="M",sigma="M",tau="M",muS="1",sigmaS="1",tauS="1",tf="1",UC="1",
     muBR="1"),
   constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001,UC=1),
   type="exgss")


p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05,
                 tf=qnorm(.1),muBR=.1,sigBR=.1,tauBR=.1)


check.p.vector(p.vector,model)
n <- c(375,375,125,125)
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,n=n,SSD=c(Inf,Inf,.25,.25)),model)

likelihood.dmc(p.vector,data)


