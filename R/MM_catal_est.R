#' Function for estimating both of the Michaelis-Menten constant and catalytic constant simultaneously
#'
#' The function estimates both of the catalytic and the Michaelis-Meten constants simultaneously using progress-curve data, 
#' enzyme concentrations, and substrate concentrations.
#' @param method This determines which model, the sQSSA or tQSSA model, is used for the estimation. Specifically, the input for method is TRUE (FALSE); then the tQSSA (sQSSA) model is used.
#' @param timespan time points when the concentrations of products were measured.
#' @param products measured concentrations of products
#' @param enz initial enzyme concentrations
#' @param subs initial substrate concentrations
#' @param K_M_init initial value of K_M constant for the Metropolis-Hastings algorithm. If the input is FALSE then it is determined by max(subs).
#' @param std standard deviation of proposal distribution. If the input is FALSE then it is determined by using the hessian of log posterior distribution.
#' @param tun tunning constant for the Metropolis-Hastings algorithm when std is FALSE (i.e., hessian of the log posterior distribution is used).
#' @param nrepeat number of effective iteration, i.e., posterior samples.
#' @param jump length of distance between sampling, i.e., thinning rate.
#' @param burn length of burn-in period.
#' @param catal_m prior mean of gamma prior for the catalytic constant k_cat.
#' @param catal_v prior variance of gamma prior for the catalytic constant k_cat.
#' @param K_M_m prior mean of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs).
#' @param K_M_v prior variance of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs)^2*1000.
#' @param volume the volume of a system. It is used to scale the product concentration. FALSE input provides automatic scaling.
#' @param t_unit the unit of time points. It can be an arbitrary string.
#' @param c_unit the unit of concentrations. It can be an arbitrary string.
#' @return A matrix containing posterior samples of the estimated parameters: the catalytic constant and the Michaelis-Menten constant.
#' @details The function MM_catal_est generates a set of Markov Chain Monte Carlo
#' simulation samples from the posterior distribution of K_M and catalytic
#' constant of enzyme kinetics model. 
#' Authors' recommendation: "Do not use this function directly. Do use the function main_est() 
#' to estimate the parameters so that the main function calls this function"
#' @examples
#' \dontrun{
#' data("timeseries_data_example")
#' timespan1=timeseries_data_example[,c(1,3,5,7)]
#' products1=timeseries_data_example[,c(2,4,6,8)]
#' MM_catal_result <- MM_catal_est(method=TRUE,timespan=timespan1,
#' products=products1,enz = c(4.4, 4.4, 440, 440), subs=c(4.4, 4.4, 4.4, 4.4), 
#' K_M_init = 1, catal_m=1, catal_v = 1000, K_M_m = 1, K_M_v = 100000, 
#' std = 10, tun =3.5, nrepeat = 1000, jump = 10, burn = 1000, 
#' volume = FALSE, t_unit = "sec", c_unit = "mM")
#' }
#' @export

MM_catal_est <- function(method, timespan, products, enz, subs, K_M_init, std, tun,
                         nrepeat, jump, burn, catal_m, catal_v, K_M_m, K_M_v, volume, t_unit, c_unit) {
  
  
  timespan = as.matrix(timespan)
  products = as.matrix(products)
  num_data = ncol(products) 
  # preprocess the imput time and product concentrations so that the product concentrations strictly increase over time.
  
  data_time <- matrix(NA, ncol = num_data, nrow = nrow(timespan))
  data_products <- matrix(NA, ncol = num_data, nrow = nrow(products))
  
  data_time[1,] = timespan[1,]
  data_products[1,] = products[1,]
  
  
  for (ii in 1:num_data){
    row_index = 1;
    for (jj in 2:sum(!is.na(timespan[,ii]))){
      if (data_products[row_index,ii] < products[jj,ii]){
        row_index = row_index + 1
        data_products[row_index,ii] = products[jj,ii]
        data_time[row_index,ii] = timespan[jj,ii]
      }
    }
  }
  max_length <- max(colSums(!is.na(data_time)))
  
  data_time = as.matrix(data_time[1:max_length,])
  data_products = as.matrix(data_products[1:max_length,])
  
  if (volume == FALSE){
    min_final <- min(apply(data_products, 2, FUN = max, na.rm =T))
    volume = (max_length - 1)/min_final
  }
  
  K_M_init = K_M_init * volume
  K_M_m = K_M_m * volume
  K_M_v = K_M_v * volume^2
  std = std * volume
  
  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal
  b_K_M = K_M_m/K_M_v
  a_K_M = K_M_m * b_K_M
  
  enz = enz * volume # scale the enzyme concentration to number.
  subs = subs * volume # scale the substrate concentration to number.
  data_products <- data_products * volume
  
  time_diff <- diff(data_time)
  prod_diff <- diff(data_products)
  S_i <- t(subs - t(data_products))
  S_i <- as.matrix(S_i[-nrow(S_i),])
  
  # log-likelihood function
  
  if (method == TRUE) {
    L.posterior <- function(km, E, S, t_diff, ni, k1, a, b) {
      num_data = ncol(S)
      n = nrow(t_diff)
      l_lik = 0
      for (jj in 1:num_data){
        for (i in 1:sum(!is.na(t_diff[,jj]))) {
          lambda_i = k1 * ((E[jj] + S[i,jj] + km)/2 - sqrt((E[jj] + S[i,jj] + km)^2 - 4 * E[jj] * S[i,jj])/2)
          l_f_i = log(dgamma(t_diff[i,jj], shape = ni[i,jj], rate = lambda_i) + 1e-300)
          l_lik = l_lik + l_f_i
        }
      }
      l_lik = l_lik + log(dgamma(km, a, b) + 1e-300)
      return(l_lik)
    }
  } else {
    L.posterior <- function(km, E, S, t_diff, ni, k1, a, b) {
      num_data = ncol(S)
      n = nrow(t_diff)
      l_lik = 0
      for (jj in 1:num_data){
        for (i in 1:sum(!is.na(t_diff[,jj]))) {
          lambda_i = (k1 * E[jj] * S[i,jj]/(km + S[i,jj]))
          l_f_i = log(dgamma(t_diff[i,jj], shape = ni[i,jj], rate = lambda_i) + 1e-300)
          l_lik = l_lik + l_f_i
        }
      }
      l_lik = l_lik + log(dgamma(km, a, b) + 1e-300)
      return(l_lik)
    }
  }

  # random walk chain sampler
  proposal <- function(km, sd) {
    ans = rnorm(1, mean = km, sd)
    while (ans < 0) {
      ans = rnorm(1, mean = km, sd)
    }
    return(ans)
  }

  MH <- function(l.star, l.m, km.star, km.m, sd) {
    logMH = l.star - l.m
    logMH = logMH + log(pnorm(km.m, mean = 0, sd)) - log(pnorm(km.star, mean = 0, sd))
    if (logMH < 0) {
      MH0 = exp(logMH)
      uni = runif(1)
      if (uni < MH0) {
        km = km.star
      } else {
        km = km.m
      }
    } else {
      km = km.star
    }
    return(km)
  }



  # Gibbs sampler of recovery parameter:

  Gibbs_catal <- function(E, S_i, km, t_diff, prod_diff, a, b) {
    num_data = ncol(S_i)
    n = nrow(t_diff)
    sum_x_t = 0
    for(jj in 1:num_data){
      if (method == T) {
        for (i in 1:sum(!is.na(t_diff[,jj]))) {
          sum_x_t = sum_x_t + ((E[jj] + S_i[i,jj] + km)/2 - sqrt((E[jj] + S_i[i,jj] + km)^2 - 4 * E[jj] * S_i[i,jj])/2) * t_diff[i,jj]
        }
      } else {
        for (i in 1:sum(!is.na(t_diff[,jj]))) {
          sum_x_t = sum_x_t + (E[jj] * S_i[i,jj]/(km + S_i[i,jj])) * t_diff[i,jj]
        }
      }
    }
    alpha_posterior = sum(prod_diff, na.rm = T) + a
    beta_posterior = sum_x_t + b
    
    k1 = rgamma(1, shape = alpha_posterior, rate = beta_posterior)
    return(k1)
  }

  ################################################################ MCMC ITERATION START
  total_repeat = burn + jump * nrepeat
  MCMC_gen = matrix(rep(0, total_repeat * 2), nrow = total_repeat, ncol = 2)
  km = K_M_init
  # catal = catal_init
  catal = Gibbs_catal(enz, S_i, km, t_diff = time_diff, prod_diff, a_catal, b_catal)
  MCMC_gen[1, ] = cbind(catal, km)
  count.km = 0
  for (rep in 2:total_repeat) {
    if (std > 0) {
      sd = std
    } else {
      hes <- hessian(func = L.posterior, km, E = enz, S = S_i, t_diff = time_diff, ni = prod_diff, k1 = catal, a = a_K_M,
                     b = b_K_M)
      sd = tun * sqrt(abs(1/hes))
      
    }

    km.star = proposal(km, sd)
    L.km.m = L.posterior(km, enz, S_i, t_diff = time_diff, prod_diff, catal, a_K_M, b_K_M)
    L.km.star = L.posterior(km.star, enz, S_i, t_diff = time_diff, prod_diff, catal, a_K_M, b_K_M)
    
    # cat(L.km.star, L.km.m, km.star, km, sd)
    
    km = MH(L.km.star, L.km.m, km.star, km, sd)
    if (km == km.star){
      count.km = count.km + 1
    }
    catal = Gibbs_catal(enz, S_i, km, t_diff = time_diff, prod_diff, a_catal, b_catal)
    MCMC_gen[rep, 1] = catal
    MCMC_gen[rep, 2] = km
    
  }

  theta = MCMC_gen[seq(from = burn+1, to = total_repeat, by = jump), ]
  theta[,2] = 1/volume * theta[,2]
  
  K_M_m = K_M_m * 1/volume
  K_M_v = K_M_v * 1/volume^2
  
  if (method == TRUE) {
    main1 = "tQSSA model: Catalytic constant (k_cat)"
    main2 = "tQSSA model: Michaelis-Menten constant (K_M)"
  } else {
    main1 = "sQSSA model: Catalytic constant (k_cat)"
    main2 = "sQSSA model: Michaelis-Menten constant (K_M)"
  }
  theta1 = theta[, 1]
  theta2 = theta[, 2]

  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta1, type = "l", main = "", xlab = "Iteration", ylab = paste("k_cat (",t_unit,"^-1)",sep = ""))
  acf(theta1, main = "")
  plot(density(theta1), main = "", xlab = paste("k_cat (",t_unit,"^-1)",sep = ""))
  mtext(side = 3, line = 1, outer = TRUE, text = main1, cex = 1.5)

  cat("MCMC simulation summary of the catalytic constant (k_cat)", "\n")
  cat("Posterior mean:       ", format(mean(theta1), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""),"\n")
  cat("Posterior sd:         ", format(sd(theta1), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""), "\n")
  cat("Credible interval(U): ", format(quantile(theta1, probs = 0.975), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""),
      "\n")
  cat("Credible interval(L): ", format(quantile(theta1, probs = 0.025), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""),
      "\n")
  cat("Relative CV:          ", format((sd(theta1)/mean(theta1))/(sqrt(catal_v)/catal_m), digits = 4, justify = "right",
                                       scientific = TRUE), "\n\n")

  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE)
  layout(mat)
  plot(theta2, type = "l", main = "", xlab = "Iteration", ylab = paste("K_M (",c_unit,")",sep = ""))
  acf(theta2, main = "")
  plot(density(theta2), main = "", xlab = paste("K_M (",c_unit,")",sep = ""))
  mtext(side = 3, line = 1, outer = TRUE, text = main2, cex = 1.5)

  cat("MCMC simulation summary of the Michaelis-Menten constant (K_M)", "\n")
  cat("Posterior mean:       ", format(mean(theta2), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""), "\n")
  cat("Posterior sd:         ", format(sd(theta2), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""), "\n")
  cat("Credible interval(U): ", format(quantile(theta2, probs = 0.975), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),
      "\n")
  cat("Credible interval(L): ", format(quantile(theta2, probs = 0.025), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),
      "\n")
  cat("Relative CV:          ", format((sd(theta2)/mean(theta2))/(sqrt(K_M_v)/K_M_m), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Acceptance ratio:     ", format(count.km/total_repeat, digits = 4, justify = "right", scientific = TRUE),
      "\n")

  par(mfrow = c(1, 1))
  if (method == TRUE) {
    plot(log10(theta[, 2]), log10(theta[, 1]), type = "p", pch = 19, cex = 0.5, main = "tQSSA model: Scatter plot", ylab = paste("Log(k_cat ( ", t_unit,"^-1))",sep = ""),
         # xlab = paste("Log(K_M (", c_unit,"))",sep = ""), xlim = c(mean(log10(theta[, 2]))-1.5, mean(log10(theta[, 2]))+1.5), ylim = c(mean(log10(theta[, 1]))-1.5, mean(log10(theta[, 1]))+1.5))
         xlab = paste("Log(K_M (", c_unit,"))",sep = ""), xlim = c(log10(mean(theta[, 2]))-1.5, log10(mean(theta[, 2]))+1.5), ylim = c(log10(mean(theta[, 1]))-1.5, log10(mean(theta[, 1]))+1.5))
  } else {
    plot(log10(theta[, 2]), log10(theta[, 1]), type = "p", pch = 19, cex = 0.5, main = "sQSSA model: Scatter plot", ylab = paste("Log(k_cat ( ", t_unit,"^-1))",sep = ""),
         xlab = paste("Log(K_M (", c_unit,"))",sep = ""), xlim = c(log10(mean(theta[, 2]))-1.5, log10(mean(theta[, 2]))+1.5), ylim = c(log10(mean(theta[, 1]))-1.5, log10(mean(theta[, 1]))+1.5))
  }
  

  theta = as.data.frame(theta)
  names(theta) = c(paste("k_cat (",t_unit,"^-1)",sep = ""), paste("K_M (",c_unit,")",sep = ""))
  # write.csv(theta, "Simul.csv")
  return(theta)
}
