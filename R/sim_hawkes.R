
#' @export
decay_func <- function(x, alpha, delta){
  if(length(x) == 0){
    return(0)
  }

  out <- alpha * delta * exp(-delta * x)
  return(out)
}

#' Simulate inhomogenous Hawkes process
#' @export
#' @param lambda_fun function for inhomogenous poisson process
#' @param history0 vector of historic times
#' @param kern function Hawkes kernal
#' @param startT numeric start time
#' @param endT numeric end time
#' @param maxI numeric maximum number of points to simulate
#' @examples
#' lambda_levels <- c(0.5, 1, 2, 1, 3, 0.5, 1, 1, 2, 0.5, 0.5, 0.5, 0.5, 2, 2, 2, 4, 0, 1, 1, 1)
#' lambda_fun <- stepfun(seq(10, 200, by = 10), y = lambda_levels)
#' kern <- function(x){return(decay_func(x, alpha = 1, delta = 1.2))}
#'
#' x <- sim_hawkes(lambda_fun, NULL, kern, 0, 200)
#'
#' hist(x, breaks = 200)
#'
#' plot(lambda_fun, add = TRUE, col = "red")
sim_hawkes <- function(lambda_fun, history0 = NULL, kern, startT = max(history0), endT = 10, maxI = Inf, progressBar = TRUE){

  TT = startT
  i = 1
  if(gtools::invalid(history0)){history0 <- NULL}
  history <- history0

  if(progressBar){
    pb <- txtProgressBar(min = startT, max = endT, style = 3)
  }

  while(i <= maxI){
    lambda_max <- lambda_fun(TT) + sum(kern(TT - history))

    u <- runif(1, 0, 1)
    tau <- -log(u) / lambda_max

    TT <- TT + tau

    s <- runif(1, 0, 1)
    lambda_current <- lambda_fun(TT) + sum(kern(TT - history))

    if(s <= lambda_current / lambda_max){
      if(TT > endT){
        return(history)
      }
      if(progressBar){setTxtProgressBar(pb, TT)}
      history <- c(history, TT)
      i <- i + 1
    }
  }

  return(history)
}
