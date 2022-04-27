
#' Generate posterior samples of static and state parameters.
#' @param prior_sample Matrix of prior samples of static parameters. The number of rows is the number of samples, the number of components of the parameter is equal to the number of columns.
#' @param dprior Function which evaluates the prior density. This function should have one argument taking the entire vector of static parameters.
#' @param loss Function which takes as input values of static and state parameters and returns a distance and an updated value for the state parameter, see Details.
#' @param loss_args List argument which is passed to the \code{loss} function. See Details.
#' @param Ntheta Positive number representing the number of static parameter proposals.
#' @param Nx Positive number representing the number of state parameter proposals for each static parameter proposal.
#' @param pacc Number between 0 and 1. The threshold for calculating the threshold for accepting state parameter proposals.
#' @param dtp Positive number. Time between time steps.
#' @param ESS_threshold Positive number. Effective sample size (ESS) threshold, when ESS drops below this number multiplied by \code{Ntheta} then static parameter proposals are replenished.
#' @param eps Positive valued numeric vector of length \code{TT}, optional. This vector is usually generated automatically from \code{pacc}, but if this argument is used then \code{pacc} is ignored.
#' @param cl An object of class "cluster" to use \code{parLapply} internally or the string "mclapply" to use \code{mclapply} internally. A NULL input (the default) will use \code{lapply} internally, a single core.
#' @param TT Positive integer. Number of time steps.
#' @param trans Function to transform parameters for resampling purposes, see Details.
#' @param invtrans Inverse function of \code{trans} to invert the transformation.
#' @param resample_times Numeric vector of replenishment times, optional.
#' @param trans_args List argument which is passed to the \code{trans} and \code{invtrans} function.
#' @param cov_coef Positive number to
#' @param acceptance_correction Positive number. Correction to acceptance probability for \code{trans} function.
#' @details
#' The \code{loss} argument is a function with the following arguments: \code{x}, \code{theta}, \code{time1}, \code{time2}, \code{inp}. The argument \code{x} is the current state parameter proposal, \code{theta} is the static parameter proposal, \code{time1} and \code{time2} are the start and end times for the update step, and \code{inp} takes a list of any other inputs to the function. This function should return list of length two with a distance and an updated value for the state parameter. See the example.
#' @examples
#' \dontrun{
#' library(ABCSMC2)
#' library(ggplot2)
#' library(ggalt)
#'
#' #cl <- makeCluster(parallel::detectCores() - 1)
#' #cl <- "mclapply"
#' cl <- NULL
#'
#' TT <- 30
#' true_theta <- c(0.25, 0.5)
#' lower <- 0
#' upper <- 3.5
#' sd_t <- 1
#' a_logit <- 0.9
#' dist_coef <- 0.5
#' true_states <-cumsum(rnorm(TT))
#'
#' lambda_fun <- stepfun(seq(1, TT - 1, by = 1), y = true_states)
#'
#' generator <- function(TT, true_states, theta){
#'     y <- list()
#'     for(i in 1:TT){
#'         obs <- as.numeric(sn::rsn(1e5, dp = sn::cp2dp(cp = c(true_states[i], theta), "SN")))
#'         y[[i]] <- c(mean(obs), stats::sd(obs), e1071::skewness(obs))
#'     }
#'     return(y)
#' }
#'
#' y <- generator(TT, true_states, true_theta)
#'
#'
#' inp <- list(
#'   lower = lower,
#'   upper = upper,
#'   sd_t = sd_t,
#'   a_logit = a_logit,
#'   y = y
#' )
#'
#' loss <- function(x, theta, time1, time2, inp){
#'     x <- as.numeric(x[1])
#'
#'     if(gtools::invalid(x)){
#'         x <- rnorm(1)
#'     } else {
#'         x <- rnorm(1, mean = x)
#'     }
#'
#'     y <- sn::rsn(1e1, dp = sn::cp2dp(cp = c(x, theta), "SN"))
#'
#'     ss_obs <- inp$y[[time2]]
#'     ss_sim <- c(mean(y), stats::sd(y), e1071::skewness(y))
#'
#'     return(list(distance = sqrt(sum((ss_obs - ss_sim)^2)), x = x))
#' }
#'
#'
#' Ntheta = 100
#' Nx = 100
#' pacc = 0.2
#'
#' lower_theta <- c(0.1, 0.2)
#' upper_theta <- c(0.5, 0.8)
#'
#' prior_sample <- data.frame(theta1 = runif(Ntheta, lower_theta[1], upper_theta[1]), theta2 = runif(Ntheta, lower_theta[2], upper_theta[2]))
#'
#' prior_sample <- as.matrix(prior_sample, ncol = 2)
#'
#' trans <- function(x, trans_args){
#'   theta1 <- qnorm(x[,1])
#'   theta2 <- qnorm(x[,2] / 0.9)
#'   return(cbind(theta1, theta2))
#' }
#'
#' invtrans <- function(x, trans_args){
#'   theta1 <- pnorm(x[,1])
#'   theta2 <- pnorm(x[,2]) * 0.9
#'   return(cbind(theta1, theta2))
#' }
#'
#' acceptance_correction <- function(x){
#'     1/dnorm(qnorm(x[1])) * (1/0.9)/(dnorm(qnorm(x[2]/0.9)))
#' }
#'
#'
#' full_list <- SMC2_ABC(prior_sample, dprior = function(x){dunif(x[1], 0.1, 0.5)*dunif(x[2],0.2,0.8)}, loss, loss_args = inp, Ntheta = Ntheta, Nx = Nx, pacc = pacc, cl = cl, dt = 1, ESS_threshold = 0.5, TT = TT, trans = trans, invtrans = invtrans, cov_coef = 0.25^2, acceptance_correction = acceptance_correction)
#'
#' state_df <- get_state(full_list)
#'
#' state_df$state <- true_states
#'
#' theta_df <- get_parameter(full_list)
#'
#'
#' ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))
#'
#' ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")
#' }
#' @export
SMC2_ABC <- function(prior_sample, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp = 1, ESS_threshold = 0.1, eps = NULL, cl = NULL, TT, trans = function(x, trans_args){I(x)}, invtrans = function(x, trans_args){I(x)}, resample_times = NA, trans_args = list(), cov_coef = 1, acceptance_correction = function(x){
  1})
{

  prior_sample <- as.matrix(prior_sample)

  # Check inputs to function
  .args <- as.list(as.list(environment()))
  do.call(check_input, .args)


  parallel <- ifelse(is.null(cl), 1,
                     ifelse("cluster" %in% class(cl), 2,
                            ifelse(cl == "mclapply", 3,
                                   ifelse(cl == "test", 4, 5)))
  )

  lfunc <- make_listfunc(parallel, cl)

  n_params <- ifelse(gtools::invalid(dim(prior_sample)[2]), 1, dim(prior_sample)[2])

  x_list <- as.list(rep(NA, Ntheta))
  x_list <- lapply(x_list, function(i){list()})

  full_list <- list()
  if(is.null(eps)){
    eps <- rep(NA, TT)
  }

  for(m in 1:Ntheta){
    x_list[[m]]$theta <- prior_sample[m, , drop = FALSE]
    x_list[[m]]$x <- matrix(NA, nrow = Nx)
    x_list[[m]]$w <- rep(1, Nx)
    x_list[[m]]$p <- NULL
    x_list[[m]]$lpprod <- NULL
    x_list[[m]]$omega <- 1
  }

  q_l <- list()


  for(tp in 1:TT){
    x_list <- lfunc(x_list, particle_filter_wrapper, time1 = tp*dtp - dtp, time2 = tp*dtp, loss = loss, loss_args = loss_args)

    distances <- sapply(x_list, function(x){x$distance})
    if(tp > 1){
      d_weights <- sapply(full_list[[tp - 1]], function(x){rep(x$omega, Nx)})
    } else {
      d_weights <- rep(1, length(distances))
    }

    if(is.na(eps[tp])){
      eps[tp] <- DescTools::Quantile(distances, d_weights, probs = pacc)
    }

    for(m in 1:Ntheta){
      x_list[[m]]$w     <- (x_list[[m]]$distance <= eps[tp])*1
      x_list[[m]]$lpprod <- sum(x_list[[m]]$lpprod, log(mean(x_list[[m]]$w)))
      x_list[[m]]$p     <- mean(x_list[[m]]$w)

      if (x_list[[m]]$p == 0) {
        x_list[[m]]$w <- rep(1, Nx)
      }
      x_list[[m]]$w <- x_list[[m]]$w / sum(x_list[[m]]$w)
      x_list[[m]]$omega <- x_list[[m]]$omega * x_list[[m]]$p
    }

    omegas <- sapply(x_list, function(x){x$omega})
    thetas <- matrix(t(sapply(x_list, function(x){x$theta})), ncol = n_params)

    ESS <- sum(omegas) ^ 2 / sum(omegas^2)
    for(i in 1:dim(thetas)[2]){
      print(as.numeric(DescTools::Quantile(thetas[,i], weights = omegas, probs = c(0, 0.025, 0.5, 0.975, 1))))
    }
    ifelse(ESS_threshold > 0, print(paste0(tp, ". ", ESS)), print(paste("SMC: (", tp, ") ", ESS)))

    full_list[[tp]] <- x_list

    if(ESS < Ntheta * ESS_threshold | tp %in% resample_times){
      print("resample")
      post_cov <- cov_coef * cov.wt(trans(thetas, trans_args = trans_args), wt=omegas * 1e20)$cov
      nb <- 0
      aa <- sample(1:Ntheta, Ntheta, prob = omegas/sum(omegas), replace = TRUE)
      proposed_log_theta <- trans(thetas[aa, , drop = FALSE], trans_args = trans_args) + mvtnorm::rmvnorm(Ntheta, sigma = post_cov)
      proposed_thetas <- invtrans(proposed_log_theta, trans_args = trans_args)

      full_list_prop <- SMC2_ABC(proposed_thetas, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp, ESS_threshold = 0, eps = eps[1:tp], cl,  TT = tp, trans = trans, invtrans = invtrans)

      x_list_prop <- full_list_prop[[tp]]

      proposed_omegas <- sapply(x_list_prop, function(x){x$omega})
      #p_aa <- sample(1:Ntheta, Ntheta, prob = proposed_omegas/sum(proposed_omegas), replace = TRUE)

      for(m in 1:Ntheta){
        #proposed_Z_hat <- x_list_prop[[p_aa[m]]]$pprod
        #old_Z_hat      <- x_list[[aa[m]]]$pprod

        Z_ratio <- exp(x_list_prop[[m]]$lpprod - x_list[[aa[m]]]$lpprod)

        MH_ratio <- Z_ratio * dprior(proposed_thetas[m,]) / (dprior(thetas[aa[m],]))

        MH_ratio <- MH_ratio * acceptance_correction(thetas[aa[m],]) / acceptance_correction(proposed_thetas[m,])

          #(0.5/(dnorm(qnorm(thetas[aa[m]])))) * (0.5/(dnorm(qnorm(proposed_thetas[p_aa[m]]))))

        un <- runif(1)

        if(un < MH_ratio){
          nb <- nb + 1
          for(time_star in 1:tp){
            full_list[[time_star]][[m]] <- full_list_prop[[time_star]][[m]]
          }
        } else {
          for(time_star in 1:tp){
            full_list[[time_star]][[m]] <- full_list[[time_star]][[aa[m]]]
          }
        }

        x_list[[m]] <- full_list[[tp]][[m]]

        x_list[[m]]$omega <- 1
      }
      cat("acceptance rate: ", nb/Ntheta, "\n")

    }


    full_list[[tp]] <- x_list

  }

  return(full_list)
}

#' @export
make_listfunc <- function(parallel, cl){
  output <- switch(parallel,
                   lapply,
                   function(X, fun, ...){
                     return(parallel::parLapply(cl, X, fun, ...))
                   },
                   parallel::mclapply
  )

  return(output)
}
