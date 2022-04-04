
#' Initialise SMC
#' @export
#' @examples
#' param <- list(alpha = 1.75, beta = 0.1, mu = -0.2, phi = 0.95, s_h = 0.6, s_v = 0.8)
#' tmp <- SMC.initialize(param, Nx, Ny, yobs[1])
#' Nx <- 1e5
Ny <- 10
SSI_initialize <- function(param, Nx, Ny, y_obs, rprior, robs) {
  tmp <- list(
    a = rep(1:Nx, times = Ny),    # not really useful here
    x = rprior(Nx, Ny),

    x = rep(rnorm(Nx), each = Ny),
    y = robs,
    dist = abs(y - y_obs)
  )
  return(tmp)
}
