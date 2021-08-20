#' Function for estimating the Michaelis-Menten constant
#'
#' The function estimates the Michaelis-Menten constant using progress-curve data, 
#' enzyme concentrations, substrate concentrations, and the catalytic constant.
#' @param method This determines which model, the sQSSA or tQSSA model, is used for the estimation. Specifically, the input for method is TRUE (FALSE); then the tQSSA (sQSSA) model is used.
#' @param timespan time points when the concentrations of products were measured.
#' @param products measured concentrations of products
#' @param enz initial enzyme concentrations
#' @param subs initial substrate concentrations
#' @param catal true value of the catalytic constant.
#' @param K_M_init initial value of K_M constant for the Metropolis-Hastings algorithm. If the input is FALSE then it is determined by max(subs).
#' @param std standard deviation of proposal distribution. If the input is FALSE then it is determined by using the hessian of log posterior distribution.
#' @param tun tuning constant for the Metropolis-Hastings algorithm when std is FALSE (i.e., hessian of the log posterior distribution is used).
#' @param nrepeat number of effective iteration, i.e., posterior samples.
#' @param jump length of distance between sampling, i.e., thinning rate.
#' @param burn length of burn-in period.
#' @param K_M_m prior mean of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs).
#' @param K_M_v prior variance of gamma prior for the Michaelis-Menten constant K_M. If the input is FALSE then it is determined by max(subs)^2*1000.
#' @param volume the volume of a system. It is used to scale the product concentration. FALSE input provides automatic scaling.
#' @param t_unit the unit of time points. It can be an arbitrary string.
#' @param c_unit the unit of concentrations. It can be an arbitrary string.
#' @return A vector containing posterior samples of the estimated parameter: the Michaelis-Menten constant.
#' @details The function MM_est generates a set of Markov Chain Monte Carlo
#' simulation samples from posterior distribution of the Michaelis-Menten constant of enzyme kinetics model. 
#' Because the function estimates only the Michaelis-Menten constant the true value of the catalytic constant should be given.
#' Authors' recommendation: "Do not use this function directly. Do use the function main_est() 
#' to estimate the parameter so that the main function calls this function"
#' @examples
#' \dontrun{
#' data("timeseries_data_example")
#' timespan1=timeseries_data_example[,c(1,3,5,7)]
#' products1=timeseries_data_example[,c(2,4,6,8)]
#' MM_result <- MM_est(method=TRUE,timespan=timespan1,products=products1,
#' enz = c(4.4, 4.4, 440, 440), subs=c(4.4, 4.4, 4.4, 4.4), catal = 0.051, 
#' K_M_init = 1, K_M_m = 1, K_M_v = 100000, std = 10, tun =3.5,
#' nrepeat = 1000, jump = 10, burn = 1000, volume = FALSE, 
#' t_unit = "sec", c_unit = "mM")
#' }
#' @export
MM_est <- function(method, timespan, products, enz, subs, catal, K_M_init, std, tun,
                   nrepeat, jump, burn, K_M_m, K_M_v, volume, t_unit, c_unit) {
    

    timespan = as.matrix(timespan)
    products = as.matrix(products)
    num_data = ncol(products)
    
    # timespan = matrix(data = c(1,1.5,2,2.3,3,3.5,3.7,4, 0.5,1,1.5,2.3,3,4,NA,NA,1.2,2,4,4.5,5,NA,NA,NA), ncol = 3)
    # products = matrix(data = c(3,3.5,5,6,7,7,8,9, 1,2,3,4,4.5,5,NA,NA,2,3,4,5,5.5,NA,NA,NA), ncol = 3) 
    # MM_est(method = T, timespan, products, enz = c(10,10,8), subs = c(30,25,25), K_M_init = 5, catal = 1, tun = 2.4, std = 10, nrepeat = 1000, jump = 10, burn = 10000)
    
    
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
    
    b_k = K_M_m/K_M_v
    a_k = K_M_m * b_k
    
    enz = enz * volume # scale the enzyme concentration to number.
    subs = subs * volume # scale the substrate concentration to number.
    data_products <- data_products * volume
    
    
    time_diff <- diff(data_time)
    prod_diff <- diff(data_products)
    S_i <- t(subs - t(data_products))
    S_i <- as.matrix(S_i[-nrow(S_i),])
    
    
    # log-likelihood function
    if (method == T) {
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


    # Trancated normal proposal
    proposal <- function(km, sd) {
        ans = rnorm(1, mean = km, sd)
        while (ans < 0) {
            ans = rnorm(1, mean = km, sd)
        }
        return(ans)
    }
    # Metropolis-Hastings step
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


    ################################################################ MCMC ITERATION START
    total_repeat = burn + jump * nrepeat
    MCMC_gen = matrix(rep(0, total_repeat), nrow = total_repeat, ncol = 1)
    
    MCMC_gen[1, ] = K_M_init
    km = K_M_init
    count.km = 0
    
    for (rep in 2:total_repeat) {
        if (std > 0) {
            sd = std
        } else {
            hes <- hessian(func = L.posterior, km, E = enz, S = S_i, t = time_diff, ni = prod_diff,
                           k1 = catal, a = a_k, b = b_k)
            sd = tun * sqrt(abs(1/hes))
        }
        km.star = proposal(km, sd)
        
        L.km.m <- L.posterior(km = km, E = enz, S = S_i, t_diff = time_diff, ni = prod_diff, k1 = catal, a = a_k, b = b_k)
        L.km.star = L.posterior(km = km.star, E = enz, S = S_i, t_diff = time_diff, ni = prod_diff, k1 = catal, a = a_k, b = b_k)
        km = MH(L.km.star, L.km.m, km.star, km, sd)
        if (km == km.star){
            count.km = count.km + 1
        }
        MCMC_gen[rep, 1] = km
    }
    
    theta = MCMC_gen[seq(from = burn+1, to = total_repeat, by = jump), ]
    theta = 1/volume * theta
    
    K_M_m = K_M_m * 1/volume
    K_M_v = K_M_v * 1/volume^2
    
    if (method == TRUE) {
        main = "tQSSA model: Michaelis-Menten constant (K_M)"
    } else {
        main = "sQSSA model: Michaelis-Menten constant (K_M)"
    }
    
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE)
    layout(mat)

    plot(theta, type = "l", main = "", xlab = "Iteration", ylab = paste("K_M (",c_unit,")",sep = ""))
    acf(theta, main = "")
    plot(density(theta), main = "", xlab = paste("K_M (",c_unit,")",sep = ""))
    mtext(side = 3, line = 1, outer = TRUE, text = main, cex = 1.5)

    cat("MCMC simulation summary of the Michaelis-Menten constant (K_M)", "\n")
    cat("Posterior mean:       ", format(mean(theta), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),"\n")
    cat("Posterior sd:         ", format(sd(theta), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),"\n")
    cat("Credible interval(U): ", format(quantile(theta, probs = 0.975), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),
        "\n")
    cat("Credible interval(L): ", format(quantile(theta, probs = 0.025), digits = 4, justify = "right", scientific = TRUE), paste("(",c_unit,")",sep = ""),
        "\n")
    cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(K_M_v)/K_M_m), digits = 4, justify = "right",
        scientific = TRUE), "\n")
    cat("Acceptance ratio:     ", format(count.km/total_repeat, digits = 4, justify = "right", scientific = TRUE),
        "\n")

    theta = as.data.frame(theta)
    names(theta) = c(paste("K_M (",c_unit,")",sep = ""))
    # write.csv(theta, "K_M.csv")
    return(theta)
}

