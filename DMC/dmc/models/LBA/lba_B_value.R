
# LSi=NULL;RSi=NULL; SRi=NULL; Ri=NULL; FBi=NULL

# LSi1=LS[i+1];RSi1=RS[i+1];LSi=LS[i]; RSi=RS[i]; Ri=out[i,"R"];
# SRi=attr(p,"SR"); FBi=cv.out[i,c("FBL","FBR")]
                     
adapt.dmc <- function(pi1,LSi1=NULL,RSi1=NULL,LSi=NULL,RSi=NULL,
  SRi=NULL,Ri=NULL,FBi=NULL) 
  # pi1 is a list of parameters for current trial
  # a is a list of quantities to be adapted (if absent no update)
  # S/FB: stimulus/response feedback (0/1) on current/previous trial
{

  # Stimulus learning function
  # SR is the representation, FB is learning signal (e.g,. feedback)
  dSR <- function(SR,FB,alpha=.1) SR + alpha*(FB-SR) 

  if ( !is.null(SRi) ) { # Update
    # Learn stimulus representation
    SRi1 <- SRi
    if ( pi1$aV[1] != 0 ) {
      if ( Ri==1 ) Ui <- LSi else Ui <- RSi
        SRi1[Ui] <- dSR(SR=SRi1[Ui],FB=as.numeric(FBi[Ri]),alpha=pi1$aV[1])
    } 
  } else SRi1 <- attr(pi1,"SR")
  # ALBA rates, difference only version 
  if ( pi1$aV[1] == 0 ) pi1$mean_v <- pi1$V0 else
    pi1$mean_v <- c(pi1$V0[1] + pi1$wV[1]*(SRi1[LSi1]-SRi1[RSi1]), 
                    pi1$V0[2] + pi1$wV[2]*(SRi1[RSi1]-SRi1[LSi1])) 
  # Calculate b
  pi1$b <- pi1$A + pi1$B 
  # Pass back updated SR
  attr(pi1,"SR") <- SRi1
  pi1
}



# save.adapt=TRUE
random.dmc <- function(p.list,model,save.adapt=FALSE)
{
  
  get.p <- function(p.list,i) {
    out <- lapply(p.list,function(x){x[,i]})
    is.SR <- unlist(lapply(strsplit(names(out),"SR"),function(x){x[1]==""}))
    attr(out,"SR") <- unlist(lapply(out[is.SR],function(x){x[1]}),use.names=FALSE)
    names(attr(out,"SR")) <- names(out[is.SR])
    out
  }
  
  n <- dim(p.list[[1]])[2]  # number of trials
  cvs <- attr(p.list,"cvs")
  facs <- attr(p.list,"facs") # factor columns + CR = numeric correct response
  LS <- as.numeric(facs$LS)
  RS <- as.numeric(facs$RS)
  # Update SR
  p <- get.p(p.list,1)
  p <- adapt.dmc(pi1=p,LSi1=LS[1],RSi1=RS[1]) # for first iteration
  out <- matrix(nrow= n, ncol=2,dimnames=list(NULL,c("RT","R")))
  if (save.adapt) {
    adapt <- vector(mode="list",length=n)
    tmp <- c(attr(p,"SR"),p[c("b","mean_v","sd_v","t0","A")]) 
    names(tmp)
    adapt[[1]] <- tmp
  }
  for (i in 1:n) {
    out[i,] <- as.numeric(rlba.norm(1,
      A=p$A,b=p$b,t0=p$t0,st0=p$st0[1],mean_v=p$mean_v,sd_v=p$sd_v,
      posdrift = attr(model,"posdrift")))

# if (out[i,"R"]==1) s=1:2 else s=2:1
# A=matrix(p$A[s],nrow=1);b=matrix(p$b[s],nrow=1)
# sd_v=matrix(p$sd_v[s],nrow=1);mean_v=matrix(p$mean_v[s],nrow=1)
# like[i] <- n1PDFfixedt0.norm(dt=out[i,"RT"]-p$t0[1],A=A,b=b,sd_v=sd_v,mean_v=mean_v,
#   posdrift=attr(model,"posdrift"))

    if ( i != n ) { # update for next iteration
      # Guide adaptation with last response
      pi1 <- get.p(p.list,i+1)
      p <- adapt.dmc(pi1,LSi1=LS[i+1],RSi1=RS[i+1],LSi=LS[i],RSi=RS[i],
                     Ri=out[i,"R"],SRi=attr(p,"SR"),FBi=cvs[i,c("FBL","FBR")])
      if (save.adapt) adapt[[i+1]] <- c(attr(p,"SR"),p[c("b","mean_v","sd_v","t0","A")]) 
    }
  }
  if (save.adapt) attr(out,"adapt") <- do.call(rbind,lapply(adapt,unlist))
  out
}


transform.dmc <- function(par.df,do.trans=TRUE) 
{
  
  if (do.trans)
    list(A=t(par.df$A),sd_v=t(par.df$sd_v),t0=t(par.df$t0),st0=t(par.df$st0),
         SR1=t(par.df$SR1),SR2=t(par.df$SR2),SR3=t(par.df$SR3),
         aV=t(par.df$aV),V0=t(par.df$V0),wV=t(par.df$wV),
         B=t(par.df$B)) else
    list(A=par.df$A,sd_v=par.df$sd_v,t0=par.df$t0,st0=par.df$st0,
         SR1=par.df$SR1,SR2=par.df$SR2,SR3=par.df$SR3,
         aV=par.df$aV,V0=par.df$V0,wV=par.df$wV,
         B=par.df$B)       
}


adapt.r.dmc <- function(pi1,Si1=NULL,a=NULL,Si=NULL,FBi=NULL) 
  # pi1 is a list of parameters for current trial
  # a is a list of quantities to be adapted 
  # S/FB: stimulus/response feedback (0/1) on current/previous trial
{

  # Stimulus learning function
  # SR is the representation, R is reponse
  dSR <- function(SR,S,alpha=.1) SR + alpha*(S-SR) 

  # Learn stimulus representation
  if (pi1[1,"aV"] != 0) {
    if (FBi==0)
      pi1[,"SR"] <-           c(dSR(SR=a[1,"SR"],S=Si,alpha=pi1[1,"aV"]),a[2,"SR"]) else
      pi1[,"SR"] <- c(a[1,"SR"],dSR(SR=a[2,"SR"],S=Si,alpha=pi1[1,"aV"]))
  }
  pi1
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10)   
{

  do.n1 <- function(x) matrix(x[attr(data,"n1.index")],ncol=dim(x)[2])
  
  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=FALSE,
      cells=attributes(data)$cells,
      cvs=data[,attr(attributes(data)$model,"cvs")],
      n1.index=attr(data,"n1.index")
  )
  
  pars <- array(unlist(p.list,use.names=FALSE),
    dim=c(dim(p.list[[1]]),length(p.list)),
    dimnames=list(NULL,NULL,names(p.list)))

  # Update SR
  for (i in 1:(dim(data)[1]-1)) {
    p <- adapt.r.dmc(pi1=pars[i+1,,c("SR","aB","aV")],
      Si1=data$stim[i+1],Si=data$stim[i],a=pars[i,,c("SR")],FBi=data[i,"FB"]-1)
    pars[i+1,,c("SR")] <- as.matrix(data.frame(p)[,c("SR")])
  }

  # Update parameters
  b <- pars[,,"A"] + pars[,,"B"]
  mean_v <- cbind(pars[,1,"V0"] + pars[,1,"wV"]*(pars[,2,"SR"]-data$stim),
                  pars[,2,"V0"] + pars[,2,"wV"]*(data$stim-pars[,1,"SR"]))


  # all(pars[,,"SR"]==adapt[,c("SR.s1","SR.s2","SR.s3")])
  # all(pars[,,"A"]==adapt[,c("A.r1","A.r2")])
  # all(pars[,,"sd_v"]==adapt[,c("sd_v.r1","sd_v.r2")])
  # all(pars[,,"t0"]==adapt[,c("t0.r1","t0.r2")])
  # all(mean_v==adapt[,c("mean_v.r1","mean_v.r2")])
  # all(b==adapt[,c("b.r1","b.r2")])
  #   
  # all(pars[,,"mean_v"]==adapt[,c("mean_v.r1","mean_v.r2")])
  # all(pars[,,"b"]==adapt[,c("b.r1","b.r2")])
  
  pmax(n1PDFfixedt0.norm(dt=data$RT-pars[,1,"t0"],
          A=do.n1(pars[,,"A"]),
          b=do.n1(b),
          mean_v=do.n1(mean_v),
          sd_v=do.n1(pars[,,"sd_v"]),
          posdrift=attr(attr(data,"model"),"posdrift")),min.like)
}



