# Template setup for simple RT LBA, B=b-A parameterization with a threhold that
# is a funciton of a covariate, ISI: b <- A + B*(1 - r*(1 - exp(-k*ISI)))
#   External parameters types: A, B, t0, mean_v, sd_v, r, k, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

transform.dmc <- function(par.df,do.trans=FALSE) 
{
  cvs <- attr(par.df,"cvs")
  if ( is.data.frame(par.df) ) {
    r <- pnorm(par.df$r[1])
    b <- par.df$A[1] + par.df$B[1]*(1 - r*(1 - exp(-par.df$k[1]*cvs$ISI))) 
    n <- dim(cvs)[1]
    n_v <- length(par.df$mean_v)
    list(b=matrix(rep(b,each=n_v),nrow=n_v),A=matrix(rep(par.df$A,times=n),nrow=n_v),
      mean_v=matrix(rep(par.df$mean_v,times=n),nrow=n_v),
      sd_v=matrix(rep(par.df$sd_v,times=n),nrow=n_v),
      t0=matrix(rep(par.df$t0,times=n),nrow=n_v),
      st0=matrix(rep(par.df$st0,times=n),nrow=n_v))
  } else {
    r <- data.frame(apply(par.df$r,2,pnorm))
    par.df$b <- par.df$A + par.df$B*(1 - r*(1-exp(-par.df$k*cvs$ISI)))
    par.df
  }
}


random.dmc <- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
              mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
              posdrift = attr(model,"posdrift"))
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10)   
{

  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=TRUE,
      cells=attributes(data)$cells,
      cvs=data[,attr(attributes(data)$model,"cvs"),drop=FALSE],
      n1.index=attr(data,"n1.index")
  )
    
  pmax(n1PDFfixedt0.norm(dt=data$RT-p.list$t0[,1],
          A=p.list$A,
          b=p.list$b,
          mean_v=p.list$mean_v,
          sd_v=p.list$sd_v,
          posdrift=attr(attr(data,"model"),"posdrift")),min.like)
}


