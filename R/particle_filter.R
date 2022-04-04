
#' @export
particle_filter <- function(theta, x, w, loss, loss_args, threshold, time1, time2, lfunc){
  a <- sample(1:I(dim(x)[1]), prob = w, replace = TRUE)
  dist_out <- apply(x[a, , drop = FALSE], MARGIN = 1, loss, theta = theta, inp = loss_args, time1 = time1, time2 = time2)
  return(dist_out)
}

