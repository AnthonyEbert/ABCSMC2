#' @export
loss_simple <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- rnorm(1)
  } else {
    x <- rnorm(1, mean = x)
  }

  y <- suppressMessages(sn::rsn(1e1, dp = sn::cp2dp(cp = c(x, theta), "SN")))

  ss_obs <- c(mean(inp$y[[time2]]), sd(inp$y[[time2]]), e1071::skewness(inp$y[[time2]]))
  ss_sim <- c(mean(y), sd(y), e1071::skewness(y))

  # ss_obs <- sn::summary(sn::selm(inp$y[[time2]]~1, opt.method = "BFGS"))@param.table[,1]
  # ss_sim <- sn::summary(sn::selm(y~1), opt.method = "BFGS")@param.table[,1]

  if(is.null(inp$weights)){
    inp$weights <- c(1,1,1)
  }

  distance_out <- sum( (ss_obs - ss_sim)^2 / inp$weights^2 )

  return(list(distance = sqrt(distance_out), x = x))
}

#' @export
generate_simple <- function(TT, true_states, theta){
  y <- list()
  for(i in 1:TT){
    y[[i]] <- as.numeric(sn::rsn(1e1, dp = sn::cp2dp(cp = c(true_states[i], theta), "SN")))
  }
  return(y)
}

#' @export
generate_stan_skew <- function(TT, true_states, theta){
  y <- list()
  for(i in 1:TT){
    y[[i]] <- as.numeric(sn::rsn(1e2, true_states[i], theta[1], theta[2]))
  }
  return(y)
}

#' @export
loss_stan_skew <- function(x, theta, time1, time2, inp){
  x <- as.numeric(x[1])

  if(gtools::invalid(x)){
    x <- rnorm(1)
  } else {
    x <- rnorm(1, mean = x)
  }

  y <- suppressMessages(sn::rsn(1e2, x, theta[1], theta[2]))

  ss_obs <- c(mean(inp$y[[time2]]), sd(inp$y[[time2]]), e1071::skewness(inp$y[[time2]]))
  ss_sim <- c(mean(y), sd(y), e1071::skewness(y))

  # ss_obs <- sn::summary(sn::selm(inp$y[[time2]]~1, opt.method = "BFGS"))@param.table[,1]
  # ss_sim <- sn::summary(sn::selm(y~1), opt.method = "BFGS")@param.table[,1]

  if(is.null(inp$weights)){
    inp$weights <- c(1,1,1)
  }

  distance_out <- sum( (ss_obs - ss_sim)^2 / inp$weights^2 )

  return(list(distance = sqrt(distance_out), x = x))
}



