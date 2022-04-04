#' @export
dHawkes <- function(theta){
  dunif(theta[1], 0.3, 0.7) * dunif(theta[2], 0.3, 0.7)
}

#' Hawkes loss function
#' @param x state numeric
#' @param theta parameter numeric
#' @param time1 first time step
#' @param time2 second time step
#' @param inp inputs
#' @examples
#' states <- generate_state(NULL, 2, 0, 3.5, sd = 1, 0.9)
#' theta <- c(0.25, 0.5)
#' TT <- length(states)
#' lambda_fun <- stepfun(seq(10, TT*10 - 10, by = 10), y = states)
#' kern <- function(x){return(decay_func(x, alpha = theta[1], delta = theta[2]))}
#' output <- sim_hawkes(lambda_fun, NULL, kern, 0, TT*10, progressBar = FALSE)
#' inp <- list(
#'   lower = 0,
#'   upper = 3.5,
#'   sd_t = 1,
#'   a_logit = 0.9,
#'   history = output,
#'   simulator = hawkes_simulator
#' )
#' loss_hawkes(states[2], theta, 10, 20, inp)
#' @export
loss_hawkes <- function(x, theta, time1, time2, inp){
  if(length(x) %in% c(0,1)){

    if(gtools::invalid(x)){
      x <- min(rgamma(1, 100, 100), inp$upper)
    } else {
      x <- generate_state(x, 1, lower = inp$lower, upper = inp$upper, sd = inp$sd, a = inp$a)
    }
  } else {
    x <- x[length(x)]
  }

  output <- dist_ss(x, theta, inp$history, time1 = time1, time2 = time2, simulator = inp$simulator)

  return(list(distance = output, x = x))
}

#' @export
dist_ss <- function(x, theta, history, time1, time2, simulator){

  if(gtools::invalid(which(history <= time1))){
    history <- NULL
  }

  sim_out <- simulator(x, theta, history[which(history <= time1 & history >= time1 - 100)], time1, time2, Ni = Inf)

  cotemp_sim <- sim_out$new_history
  cotemp_obs <- history[which(history > time1 & history < time2)]

  sim_ss <- sum_stat(cotemp_sim, time1, time2)
  obs_ss <- sum_stat(cotemp_obs, time1, time2)

  # est1 <- c(0, 1, -4.12e-01, 1.39e-01, 0        , 0        , 0       ,  5.29e-04, -3.17e-05, 0)
  # est2 <- c(0, 1, 0        , 1.72e-01,  1.80e-01, 0        , 0       , -3.33e-04,  2.84e-05, 0)
  # est3 <- c(0, 1, -5.34e-01, 3.23e-01, -1.07e-01, -5.08e-02, 6.09e-02,         0,         0, 0)

  est1 <- c(1.388739e-01, 4.997625e-03, 1.209265e-06, -4.838022e-07, 3.766624e-09, -1.145458e-11, 1.249928e-14, 2.635814e-03, -1.577590e-04, -1.021220e-03)
  est2 <- c(6.417910e-01, 2.508705e-03, -2.210345e-04, 9.340399e-07, 2.108424e-08, -1.946335e-10, 4.520269e-13, 2.967803e-03, -2.532569e-04, 2.830355e-03)
  est3 <- c(1.802715e-01, 5.490010e-02, 2.827284e-04, -2.936200e-05, 3.912895e-07, -2.008020e-09, 3.616250e-12,  4.192238e-04, -5.677054e-05, 5.512495e-03)

  est_mat <- rbind(est1, est3)

  sim_est <- est_mat %*% sim_ss
  obs_est <- est_mat %*% obs_ss

  return(sqrt(sum((sim_est - obs_est)^2)))
}

#' @export
dist_h <- function(x, theta, history, time1, time2, simulator){

  if(gtools::invalid(which(history <= time1))){
    history <- NULL
  }

  sim_out <- simulator(x, theta, history[which(history <= time1 & history >= time1 - 100)], time1, time2, Ni = Inf)

  cotemp_history <- history[which(history > time1 & history < time2)]

  n_out <- (sim_out$n - length(cotemp_history))^2 / max(1, length(cotemp_history))

  # true_input <- sort(c(cotemp_history, time2, time1))
  # sim_input <- sort(c(sim_out$new_history, time2, time1))
  #
  # diff_out <- mean(diff(true_input[which(!is.na(true_input))])) - mean(diff(sim_input[which(!is.na(sim_input))]))

  # if(is.nan(diff_out)){
  #   kol <- 2+2
  # }
  #
  dist_out <- n_out


  return(dist_out)
}

#' @export
sum_stat <- function(events, time1, time2){
  output <- as.numeric(na.exclude(events))
  n_events <- length(output)

  output <- c(time1, output, time2)
  diffs2    <- sum((diff(output))^2)
  diffs3     <- sum((diff(output))^3)
  min_diff <- min(diff(output))

  output <- c(int = 1, n = n_events, n2 = n_events^2, n3 = n_events^3, n4 = n_events^4, n5 = n_events^5, n6 = n_events^6, diffs2 = diffs2, diffs3 = diffs3, min_diff = min_diff)

  if(anyNA(output)){
    print(output)
  }

  return(c(int = 1, n = n_events, n2 = n_events^2, n3 = n_events^3, n4 = n_events^4, n5 = n_events^5, n6 = n_events^6, diffs2 = diffs2, diffs3 = diffs3, min_diff = min_diff))
}


#' Generate the state
#' @param init initial value
#' @param n number of states to generate
#' @param lower lower bound
#' @param upper upper bound
#' @param sd standard deviation of normal
#' @param a coefficient to apply
#' @export
#' @examples
#' init <- 2.5
#' x <- generate_state(init, 10, 0, 3.5, 1, 0.9)
generate_state <- function(init, n, lower, upper, sd, a = 1){

  if(gtools::invalid(init)){
    init <- rgamma(1, 10, 10)
    if(n == 1){
      return(init)
    }
  }

  one_step <- FALSE

  if(n == 1){
    one_step <- TRUE
    n <- 2
  }

  x <- rep(NA, n)
  x[1] <- init

  for(t in 2:n){
    x[t] <- gtools::inv.logit(a * gtools::logit(x[t-1], min = lower, max = upper) + rnorm(1, sd = sd), min = lower, max = upper)
  }

  if(one_step){
    x <- x[2]
  }

  return(x)
}
