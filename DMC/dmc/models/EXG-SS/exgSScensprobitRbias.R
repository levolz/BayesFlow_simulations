# stop-signal context independent (seperate go and stop) parameterization with tf and
# gf on the probit scale
#    External parameters types: qnorm(tf), qnorm(gf), mu, sigma, tau, muS, sigmaS, tauS 
#    Internal parameters types: tf, gf, mu, sigma, tau 

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

my.integrate <- function(...,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (class(out)=="try-error") 0 else 
  {
    if (out$subdivisions==1) 
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (class(out)=="try-error") 0 else
      {
        if (out$subdivisions==1) 0 else out$value   
      }
    } else out$value
  }
}

transform.dmc <- function(par.df) 
  # This function transfroms parameters to a form suitbale for the model 
  #   being used. Called inside of get.par.mat. 
  # "par.df" is a data frame of parameters types , some of which may need to be 
  #   transformed, or new columns created, so that the full set of internal 
  #   parameter types, specified in "type.par.names", required by the type of 
  #   evidence accumulation model being used is present.
{
  if (dim(par.df)[1]>0) {
    if (any(names(par.df)=="muBR")) par.df["RIGHT","mu"] <- par.df["RIGHT","mu"]*par.df[1,"muBR"]
    if (any(names(par.df)=="sigBR")) par.df["RIGHT","sigma"] <- par.df["RIGHT","sigma"]*par.df[1,"sigBR"]
    if (any(names(par.df)=="tauBR")) par.df["RIGHT","tau"] <- par.df["RIGHT","tau"]*par.df[1,"tauBR"]
  }
  # copy stop parameters onto go for stop accumulator
  par.df["NR",c("mu","sigma","tau")] <- par.df["NR",c("muS","sigmaS","tauS")]
  par.df$tf <- pnorm(par.df$tf)
  par.df[,c("mu","sigma","tau","muS","sigmaS","tauS","tf","UC")]
}

random.dmc<- function(n,p.df,model,SSD=Inf,staircase=NA,TRIALS=NULL)
{
  rexgss.cens(n,mu=p.df$mu,sigma=p.df$sigma,tau=p.df$tau,
         tf=p.df$tf[1],
         SSD=SSD,staircase=staircase,UpperCensor=p.df$UC[1])  
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
  # Returns vector of likelihoods for each RT in data (in same order)
{
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF.exgss.cens(dt=data$RT[attr(data,"cell.index")[[i]]],
                  tau=p.df$tau, 
                  mu=p.df$mu,
                  sigma=p.df$sigma,
                  # Trigger failure
                  tf=p.df$tf[1],
                  UpperCensor=p.df$UC[1],
                  # Stop-signal delays
                  SSD=data$SSD[attr(data,"cell.index")[[i]]],
                  # Index of stop signal accumulator
                  Si=c(1:dim(p.df)[1])[row.names(p.df)=="NR"]
      )
  }
  pmax(likelihood,min.like)
}


