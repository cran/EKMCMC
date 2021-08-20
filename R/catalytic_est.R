#' Function for estimating the catalytic constant
#'
#' The function estimates catalytic constant using progress-curve data, 
#' enzyme concentrations, substrate concentrations, and the Michaelis-Meten constant.
#' @param method This determines which model, the sQSSA or tQSSA model, is used for the estimation. Specifically, the input for method is TRUE (FALSE); then the tQSSA (sQSSA) model is used.
#' @param timespan time points when the concentrations of products were measured.
#' @param products measured concentrations of products
#' @param enz initial enzyme concentrations
#' @param subs initial substrate concentrations
#' @param K_M true value of the Michaelis-Menten constant.
#' @param catal_m prior mean of gamma prior for the catalytic constant k_cat.
#' @param catal_v prior variance of gamma prior for the catalytic constant k_cat.
#' @param nrepeat number of effective iteration, i.e., posterior samples.
#' @param jump length of distance between sampling, i.e., thinning rate.
#' @param burn length of burn-in period.
#' @param volume the volume of a system. It is used to scale the product concentration. FALSE input provides automatic scaling.
#' @param t_unit the unit of time points. It can be an arbitrary string.
#' @param c_unit the unit of concentrations. It can be an arbitrary string.
#' @return A vector containing posterior samples of the estimated parameter: the catalytic constant.
#' @details The function catalytic_est generates a set of Monte Carlo
#' simulation samples from posterior distribution of the catalytic constant of enzyme kinetics model. 
#' Because the function estimates only the catalytic constant, the true value of the Michaelis-Menten constant should be given.
#' Authors' recommendation: "Do not use this function directly. Do use the function main_est() 
#' to estimate the parameter so that the main function calls this function"
#' 
#' @examples
#' \dontrun{
#' data("timeseries_data_example")
#' timespan1=timeseries_data_example[,c(1,3,5,7)]
#' products1=timeseries_data_example[,c(2,4,6,8)]
#' catalytic_result <- catalytic_est(method=TRUE,timespan=timespan1,
#' products=products1,enz = c(4.4, 4.4, 440, 440), subs=c(4.4, 4.4, 4.4, 4.4), 
#' K_M=44, catal_m = 1, catal_v = 1000, jump = 10, burn = 1000, nrepeat = 1000, 
#' volume = FALSE, t_unit = "sec", c_unit = "mM")
#' }
#' @export
catalytic_est <- function(method, timespan, products, enz, subs, K_M,
                          catal_m, catal_v, nrepeat, jump, burn, volume, t_unit, c_unit){
 
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
    
    # dat = cbind(time, products * volume)
    K_M = K_M * volume
    enz = enz * volume # scale the enzyme concentration to number.
    subs = subs * volume # scale the substrate concentration to number.
    data_products <- data_products * volume
    
    b_k = catal_m/catal_v
    a_k = catal_m * b_k
    
    
    time_diff <- diff(data_time)
    prod_diff <- diff(data_products)
    S_i <- t(subs - t(data_products))
    S_i <- as.matrix(S_i[-nrow(S_i),])
    
    ################################################################ MCMC ITERATION START
    total_repeat = burn + jump * nrepeat
    MCMC_gen = matrix(rep(0, total_repeat), nrow = total_repeat, ncol = 1)
    
    # if (catal_init == FALSE){
    #   catal_init = 1 # this should be adjusted to a more relevant value.
    # }
    # 
    # MCMC_gen[1, ] = catal_init
    # 
    # catal = catal_init
    # 
    sum_x_t = 0
    
    for(jj in 1:num_data){
        if (method == TRUE) {
            for (i in 1:sum(!is.na(time_diff[,jj]))) {
                sum_x_t = sum_x_t + ((enz[jj] + S_i[i,jj] + K_M)/2 - sqrt((enz[jj] + S_i[i,jj] + K_M)^2 - 4 * enz[jj] * S_i[i,jj])/2) * time_diff[i,jj]
            }
        } else {
            for (i in 1:sum(!is.na(time_diff[,jj]))) {
                sum_x_t = sum_x_t + (enz[jj] * S_i[i,jj]/(K_M + S_i[i,jj])) * time_diff[i,jj]
            }
        }
    }
    
    
    alpha_posterior = sum(prod_diff, na.rm = TRUE) + a_k
    beta_posterior = sum_x_t + b_k
    
    for (rep in 1:total_repeat) {
        catal = rgamma(1, shape = alpha_posterior, rate = beta_posterior)
        MCMC_gen[rep, 1] = catal
    }
    
    # We know the exact posterior distribution of the catalytic constant.
    # In fact, we don't need to sample from the distribution. We can present the mean and the variance of the posterior distribution.
    
    theta = MCMC_gen[seq(from = burn+1, to = total_repeat, by = jump), ]

    if (method == TRUE) {
        main = "tQSSA model: Catalytic constant (k_cat)"
    } else {
        main = "sQSSA model: Catalytic constant (k_cat)"
    }
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE)
    layout(mat)
    plot(theta, type = "l", main = "", xlab = "Iteration", ylab = paste("k_cat (",t_unit,"^-1)",sep = ""))
    acf(theta, main = "")
    plot(density(theta), main = "", xlab = paste("k_cat (",t_unit,"^-1)",sep = ""))
    mtext(side = 3, line = 1, outer = TRUE, text = main, cex = 1.5)

    cat("MCMC simulation summary of the catalytic constant", "\n")
    cat("Posterior mean:       ", format(mean(theta), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""), "\n")
    cat("Posterior sd:         ", format(sd(theta), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""), "\n")
    cat("Credible interval(U): ", format(quantile(theta, probs = 0.975), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""),
        "\n")
    cat("Credible interval(L): ", format(quantile(theta, probs = 0.025), digits = 4, justify = "right", scientific = TRUE), paste("(",t_unit,"^-1)",sep = ""),
        "\n")
    cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(catal_v)/catal_m), digits = 4, justify = "right",
        scientific = TRUE), "\n")

    theta = as.data.frame(theta)
    names(theta) = c(paste("k_cat (",t_unit,"^-1)",sep = ""))
    # write.csv(theta, "catalytic.csv")
    return(theta)
}

