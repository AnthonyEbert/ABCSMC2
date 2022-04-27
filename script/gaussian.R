

# Run gaussian example which doesn't need a correction, and one that does
# The results should be exactly the same

library(ABCSMC2)
library(ggplot2)
library(ggalt)
library(dplyr)

sessionInfo()


set.seed(1)

cl <- NULL

TT <- 30
true_theta <- log(2) # logarithm of sigma

true_theta

true_states <-cumsum(rnorm(TT))

generator <- function(TT, true_states, theta){
  y <- list()
  for(i in 1:TT){
    obs <- rnorm(1e2, mean = true_states[i], sd = exp(theta[1]))
    y[[i]] <- c(mean(obs), stats::sd(obs))
  }
  return(y)
}

y <- generator(TT, true_states, true_theta)

inp <- list(
  y = y
)

loss <- function(x, theta, time1, time2, inp){

  if(!is.numeric(x)){
    x <- rnorm(1)
  } else {
    x <- rnorm(1, mean = x)
  }

  y <- rnorm(1e2, mean = x, sd = exp(theta[1]))

  ss_obs <- inp$y[[time2]]
  ss_sim <- c(mean(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE))

  return(list(distance = as.numeric(sqrt(c(1, 1) %*% (ss_obs - ss_sim)^2)), x = x))
}


Ntheta = 100
Nx = 100
pacc = 0.5

# Log sigma distribution
# Gaussian example where no correction is needed in sampler

prior_sample_log <- data.frame(theta1 = rnorm(Ntheta, 0, 5))

prior_sample_log <- as.matrix(prior_sample_log, ncol = 1)

set.seed(2)

full_list_log <-
  SMC2_ABC(
    prior_sample_log,
    dprior = function(x) {
      dnorm(x[1],0, 5)
    },
    loss,
    loss_args = inp,
    Ntheta = Ntheta,
    Nx = Nx,
    pacc = pacc,
    cl = cl,
    dt = 1,
    ESS_threshold = 0.2,
    TT = TT,
    cov_coef = 4
  )

state_df <- get_state(full_list_log)

state_df$state <- true_states

theta_log_df <- get_parameter(full_list_log)


ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_log_df[which(theta_log_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")


# Sigma distribution
# Log Gaussian example where correction is needed in sampler

loss <- function(x, theta, time1, time2, inp){

  if(!is.numeric(x)){
    x <- rnorm(1)
  } else {
    x <- rnorm(1, mean = x)
  }

  y <- rnorm(1e2, mean = x, sd = theta[1])

  ss_obs <- inp$y[[time2]]
  ss_sim <- c(mean(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE))

  return(list(distance = as.numeric(sqrt(c(1, 1) %*% (ss_obs - ss_sim)^2)), x = x))
}




prior_sample <- data.frame(theta1 = rlnorm(Ntheta, 0, 5))

prior_sample <- as.matrix(prior_sample, ncol = 1)

prior_sample <- exp(prior_sample_log)

trans <- function(x, trans_args){
  theta1 <- log(x)
  return(theta1)
}

invtrans <- function(x, trans_args){
  theta1 <- exp(x)
  return(theta1)
}

acceptance_correction <- function(x){
  1/x
}

set.seed(2)

full_list <-
  SMC2_ABC(
    prior_sample,
    dprior = function(x) {
      dlnorm(x[1],0, 5)
    },
    loss,
    loss_args = inp,
    Ntheta = Ntheta,
    Nx = Nx,
    pacc = pacc,
    cl = cl,
    dt = 1,
    ESS_threshold = 0.2,
    TT = TT,
    cov_coef = 4,
    trans = trans,
    invtrans = invtrans,
    acceptance_correction = acceptance_correction
  )

state_df <- get_state(full_list)

state_df$state <- true_states

theta_df <- get_parameter(full_list)


ggplot(state_df) + aes(x = time, y = state, ymin = lower, ymax = upper) + geom_step() + geom_ribbon(alpha = 0.2, stat = "stepribbon", fill = "red") + geom_step(mapping = aes(x = time, y = med), col = "red") + ggthemes::theme_base() + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0))

ggplot(theta_df[which(theta_df$time %% 5 == 0),]) + aes(x = value, weights = weight, col = factor(time)) + geom_density() + facet_wrap(~parameter, scales = "free")

theta_df_combined <- dplyr::bind_rows(theta_log_df %>% filter(time == TT) %>% mutate(value = exp(value)), theta_df %>% filter(time == TT), .id = "groups")

ggplot(theta_df_combined) + aes(x = value, weights = weight) + geom_density() + facet_wrap(~parameter + groups, scales = "free")
