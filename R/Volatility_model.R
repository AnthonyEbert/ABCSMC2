# Volatility model

#' @export
loss_volatility <- function(x, theta, time1, time2, inp){
  phi <- theta[1]

  if(gtools::invalid(x)){
    x <- rnorm(1, inp$mu/(1-phi), sd = sqrt(inp$s_h^2/(1-phi^2)))
  } else {
    x <- inp$mu + phi * x + inp$s_h * rnorm(1)
  }

  if(is.null(inp$gamma)){
    inp$gamma <- 0
  }

  y_prop <- exp(x/2) * stabledist::rstable(1, inp$alpha, inp$beta, inp$gamma, inp$s_v)
  output <- abs(inp$y[time2] - y_prop)

  return(list(distance = output, x = x))
}

#' @export
rprior_vol <- function(N) {
  phi <- runif(N, min = -1, max = +1)
  return(phi)
}

#' @export
dprior_vol <- function(theta){
  dunif(theta[1], -1, 1)
}
