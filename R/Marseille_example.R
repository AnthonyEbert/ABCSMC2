

#' @export
data_simulator <- function(theta, TT){

  x <- rep(NA, TT)
  y <- rep(NA, TT)

  x[1] <- rnorm(1, 1, 2)

  for(t in 2:TT){
    x[t] <- x[t-1] + rnorm(1, 0, theta[1])
  }

  y <- rnorm(TT, x, theta[2])

  return(y)
}

#' @export
simple_simulator <- function(x, theta){

  y <- x + rnorm(length(x), 0, theta)

  return(y)
}

#' @export
hawkes_simulator <- function(state, theta, history0, time1, time2, Ni = 1000){

  kern <- function(x){return(decay_func(x, alpha = theta[1], delta = theta[2]))}

  lambda_fun <- function(t){return(1*state)}

  history <- sim_hawkes(lambda_fun, history0 = history0, kern = kern, startT = time1, endT = time2, progressBar = FALSE, maxI = Ni)

  if(is.null(history)) history <- NA

  n <- length(history) - length(history0)

  new_history <- setdiff(history, history0)

  out <- list(history = history, n = n, new_history = new_history)

  return(out)
}









