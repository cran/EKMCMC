#' Main function for estimating catalytic constant k_cat and Michaelis-Menten (MM) constant K_M
#'
#' The function estimates either the catalytic constant, the Michaelis-Menten constant, or both simultaneously using progress-curve data, 
#' initial enzyme concentrations, and initial substrate concentrations.
#' @param method This determines which model, the sQSSA or tQSSA model, is used for the estimation. Specifically, the input for method is TRUE (FALSE); then the tQSSA (sQSSA) model is used. Its default value is TRUE.
#' @param timeseries Data frame containing the time points and measured concentrations of products. Every two columns represent the time points when the concentrations of the products were measured and the corresponding measured concentrations. 
#' @param enz initial enzyme concentrations
#' @param subs initial substrate concentrations
#' @param K_M true value of the Michaelis-Menten constant. Specify this object if the true value is known. Its default value is FALSE.
#' @param catal true value of the catalytic constant.  Specify this object if the true value is known. Its default value is FALSE.
#' @param K_M_init initial value of K_M constant for the Metropolis-Hastings algorithm. If the input is FALSE then it is determined by max(subs). Its default value is FALSE.
#' @param std standard deviation of proposal distribution. If the input is FALSE then it is determined by using the hessian of log posterior distribution. Its default value is FALSE.
#' @param tun tuning constant for the Metropolis-Hastings algorithm when std is FALSE (i.e., hessian of the log posterior distribution is used). Its default value is 2.4.
#' @param nrepeat number of effective iteration, i.e., posterior samples. Its default value is 1,000.
#' @param jump length of distance between sampling, i.e., thinning rate. Its default value is 10.
#' @param burn length of burn-in period. Its default value is 1,000.
#' @param catal_m prior mean of gamma prior for the catalytic constant k_cat. Its default value is 1.
#' @param catal_v prior variance of gamma prior for the catalytic constant k_cat Its default value is 1e+06.
#' @param K_M_m prior mean of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs). Its default value is FALSE.
#' @param K_M_v prior variance of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs)^2*1000. Its default value is FALSE.
#' @param volume the volume of a system. It is used to scale the product concentration. FALSE input provides automatic scaling. Its default value is FALSE.
#' @param t_unit the unit of time points. It can be an arbitrary string.
#' @param c_unit the unit of concentrations. It can be an arbitrary string.
#' @return A vector (or matrix) containing posterior samples of the estimated parameter(s).
#' @details The function main_est generates a set of Markov Chain Monte Carlo (MCMC)
#' simulation samples from the posterior distribution of the catalytic constant or (and)
#' the Michaelis-Menten constant of enzyme kinetics model. Users should
#' input initial enzyme concentrations, substrate concentrations, and progress-curve data.
#' Prior information for both parameters can be given.
#' The Gibbs sampling and Metropolis Hastings algorithms are used to sample the parameters. 
#' Parameters for the MCMC such as tuning parameter for proposal distribution, prior parameters, 
#' and the iteration number can be specified by users. This function use one of 
#' catalytic_est(), MM_est(), MM_catal_est() to generate the samples depending on parameter(s) to be estimated.
#' @importFrom numDeriv hessian
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf density quantile sd rexp
#' @importFrom graphics par layout mtext plot lines
#' @examples
#' \dontrun{
#' data("timeseries_data_example")
#' result <- main_est(method=TRUE, timeseries = timeseries_data_example, 
#' enz = c(4.4, 4.4, 440, 440), subs=c(4.4, 4.4, 4.4, 4.4), K_M_init = 1e+1, 
#' std=1e+1, tun = 3.5, jump=10, burn=1000, nrepeat=1000,
#' catal_m=1, catal_v=100, K_M_m=1, K_M_v=1e+4, volume = FALSE, 
#' t_unit = "sec", c_unit = "mM")
#' }
#' @export

main_est <- function(method = TRUE, timeseries, enz, subs, K_M = FALSE, catal = FALSE, K_M_init = FALSE, std = FALSE, tun = 2.4,
                     nrepeat = 1000, jump = 10, burn = 1000, catal_m = 1, catal_v = 1e+06, K_M_m = FALSE, K_M_v = FALSE, volume = FALSE, t_unit, c_unit) {
  num_p_curves <- dim(timeseries)[2]/2;
  timespan <- timeseries[, seq(from = 1, to=2*num_p_curves, by = 2)];
  products <- timeseries[, seq(from = 2, to=2*num_p_curves, by = 2)];
  enz <- as.numeric(enz)
  subs <- as.numeric(subs)
  
  if (K_M_init == FALSE){
    K_M_init = max(subs)
  }
  if (K_M_m == FALSE){
    K_M_m = max(subs)
  }
  if (K_M_v == FALSE){
    K_M_v = max(subs)^2 * 1000
  }
  
  if (K_M > 0 & catal == FALSE){
    result1 <- catalytic_est(method = method, timespan = timespan, products = products,enz = enz, subs = subs, 
                  K_M = K_M, nrepeat = nrepeat, jump = jump, burn = burn, catal_m=catal_m, catal_v=catal_v,
                  volume = volume, t_unit = t_unit, c_unit = c_unit)
  }else if (K_M == FALSE & catal > 0){
    result1 <-MM_est(method = method, timespan = timespan, products = products, enz = enz, subs = subs,
           K_M_init = K_M_init, catal = catal, tun = tun, std = std,
           nrepeat = nrepeat, jump=jump,burn=burn,K_M_m=K_M_m,K_M_v=K_M_v,volume=volume, t_unit = t_unit, c_unit = c_unit)
  }else if (K_M == FALSE & catal == FALSE){
    result1 <- MM_catal_est(method=method, timespan = timespan,products = products,enz=enz,subs=subs,
                 K_M_init=K_M_init,tun=tun,std=std,nrepeat=nrepeat,jump=jump,burn=burn,
                 catal_m=catal_m,catal_v=catal_v, K_M_m=K_M_m,K_M_v=K_M_v,volume=volume, t_unit = t_unit, c_unit = c_unit)
  }else{
    warning("Input true parameters values or initial values are inappropriate.\n")
    reset1 <- FALSE
  }
  return(result1)
}




