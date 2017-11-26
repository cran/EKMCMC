#' Simulation plot of enzyme kinetics model
#'
#' The function depicts the overlayed two plots; one is observed data,
#' the other is simulation result using fitted MM constant and catalytic constant.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param CL Adding empircal 95\% confidence interval (default = T)
#' @param time observed time interval
#' @param species observed trajectory of product
#' @param enz enzyme concentration
#' @param subs substrate concentration
#' @param MM true value of MM constant
#' @param catal initial value of catalytic constant
#' @param nrepeat total number of simulation (default=100)
#' @param ti tme interval for descreted simulatin result (default =1)
#' @return This functio has no returned object.
#' @details Basically this function draws overlayed picture: The one is trajectory
#' of given data of products for enzyme kinetics model. The other is trajectory of
#' products from simulation result of Gillespie algorithm with estimated two
#' constants, Michaelis-Menten constant and catalytic constant.
#' CL option controls the plot that depicts the observed data and
#' mean of simulated series only or adding 95% empirical
#' confidence interval with 10 samples of simulated trajectory.
#' @examples
#' data("Chymo_low")
#' time1=Chymo_low[,1]
#' species1=Chymo_low[,2]
#' fore_p(method=TRUE, CL=TRUE, time=time1,species=species1,enz=4.4e+7, subs=4.4e+7
#'        ,MM=4.4e+8, catal=.051)
#' @export
fore_p <- function(method = TRUE, CL=TRUE, time, species, enz, subs, MM, catal
                   , nrepeat = 100, ti = 1){
  scale = (length(time) - 1)/subs
  dat = cbind(time, species * scale)
  enz = enz * scale
  subs = subs * scale
  MM = MM * scale

  if(method==TRUE){
    Gillespie <- function(E=enz, S=subs, kd=MM, k1=catal){
      tvec=0;t=0;
      while (S>0){
        lambda = k1 * (E + kd + S - sqrt((E + kd +S)^2 - 4*E*S))/2
        t = t + rexp(1,rate=lambda)
        tvec = rbind(tvec, t)
        S = S - 1
      }
      xmat = 0:(length(tvec)-1)
      return(list(t=tvec, P=xmat))
    }
  }else{
    Gillespie <- function(E=enz, S=subs, kd=MM, k1=catal){
      tvec=0;t=0;
      while (S>0){
        lambda = (k1 * E * S)/(kd + S)
        t = t + rexp(1,rate=lambda)
        tvec = rbind(tvec, t)
        S = S - 1
      }
      xmat = 0:(length(tvec)-1)
      return(list(t=tvec, P=xmat))
    }
  }

  discretise <-function(out, dt=1){
    start = 0
    events = length(out$t)
    end = out$t[events]
    len = (end-start)%/%dt+1
    x = matrix(0,nrow=len,ncol=1)
    tvec = matrix(0,nrow=len,ncol=1)
    target = 0
    j=1
    for (i in 1:events){
      while (out$t[i]>=target){
        x[j]=out$P[i]
        tvec[j]=target
        j=j+1
        target=target+dt
      }
    }
    return(cbind(tvec,x))
  }

  if(CL == TRUE){
    par(mfrow=c(1,1))
    plot(time, species, pch=18, xlab="Time", ylab="Product")
    out1 <- Gillespie(enz,subs,MM,catal)
    fit <- discretise(out1, dt=ti)
    colnames(fit) <- c("tvec", " ")
    for(i in 2:nrepeat){
      out1 <- Gillespie(enz,subs,MM,catal)
      fit1 <- discretise(out1, dt=ti)
      colnames(fit1) <- c("tvec", " ")
      fit <- merge(fit, fit1, by ="tvec", all=T )
      if(i%%10==0) lines(fit1[,1],fit1[,2]/scale, col="grey")
    }
    fit.t1 <- fit[,1]
    fit.P <- fit[,2:(nrepeat +1)]
    fit.P1 <- apply(fit.P, 1, function(x) mean(x,na.rm = TRUE))
    fit.l  <- apply(fit.P, 1, function(x) quantile(x,0.025,na.rm = TRUE))
    fit.u  <- apply(fit.P, 1, function(x) quantile(x,0.975,na.rm = TRUE))
    lines(fit.t1, fit.l/scale,ty='l',lty=1, col='blue')
    lines(fit.t1, fit.u/scale,ty='l',lty=1, col='blue')
    lines(fit.t1, fit.P1/scale, pch=18, col="red")

  }else{
    out1 <- Gillespie(enz,subs,MM,catal)
    fit <- discretise(out1, dt=ti)
    colnames(fit) <- c("tvec", " ")
    for(i in 2:nrepeat){
      out1 <- Gillespie(enz,subs,MM,catal)
      fit1 <- discretise(out1, dt=ti)
      colnames(fit1) <- c("tvec", " ")
      fit <- merge(fit, fit1, by ="tvec", all=T )
    }
    fit.t1 <- fit[,1]
    fit.P <- fit[,2:(nrepeat+1)]
    fit.P1 <- apply(fit.P, 1, function(x) mean(x,na.rm = TRUE))
    par(mfrow=c(1,1))
    plot(time, species, pch=18, xlab="Time", ylab="Product")
    lines(fit.t1, fit.P1/scale, pch=18, col="red")
  }
}




