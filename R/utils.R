
#' Check input to SMC2_ABC function
#' @noRd
check_input <- function(prior_sample, dprior, loss, loss_args, Ntheta, Nx, pacc, dtp = 1, ESS_threshold = 0.1, eps = NULL, cl = NULL, TT, trans = I, invtrans = I, resample_times = NA, trans_args = list(), cov_coef = 1, acceptance_correction = function(x){
  1/(prod(x))
}){

  stopifnot(as.matrix(is.matrix(prior_sample)))
  stopifnot(is.function(dprior))

  stopifnot(length(formals(loss)) == 5)
  stopifnot(all(c("x", "theta", "time1", "time2", "inp") %in% names(formals(loss))))

  loss_out <- loss(NULL, prior_sample[1,], 1 * dtp, 2 * dtp, loss_args)

  stopifnot(is.list(loss_out))

  #stopifnot(is.list(loss_args))
  stopifnot(Ntheta > 0 & Ntheta %% 1 == 0)
  stopifnot(Nx > 0 & Ntheta %% 1 == 0)
  stopifnot(pacc >= 0 & pacc <= 1)
  stopifnot(dtp > 0)
  stopifnot(is.numeric(ESS_threshold))
  stopifnot(is.null(eps) | length(eps) == TT)

  stopifnot(is.function(trans))
  stopifnot(is.function(invtrans))

  # stopifnot(length(formals(trans)) == 2)
  stopifnot(length(formals(trans)) == length(formals(invtrans)))

  if(
    isFALSE(all.equal(
      prior_sample[1, , drop = FALSE],
      invtrans(
        trans(prior_sample[1, , drop = FALSE], trans_args),
      trans_args)
    )
  )) {

    warning("invtrans does not undo trans function")
  }

  stopifnot(cov_coef > 0)
  stopifnot(is.function(acceptance_correction))


  return(TRUE)

}
