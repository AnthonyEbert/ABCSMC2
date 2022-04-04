
#' @export
theta_filter <- function(input_list, time1, time2, loss = loss, loss_args){

  output <- particle_filter(input_list$theta, input_list$x, input_list$w, loss = loss, loss_args = loss_args, time1 = time1, time2 = time2)

  size_x <- length(output[[1]]$x)

  output_list <- input_list

  output_list$distance <- sapply(output, function(x){x$distance})

  output_list$x <- matrix(sapply(output, function(x){x$x}), byrow = TRUE, ncol = size_x)

  return(output_list)
}

#' @export
weighter_filter <- function(input_list, eps){

  output <- lapply(input_list$distances)

  distances <- input_list$distance

  output <- (distances <= eps) * 1

  input_list$w <- output

  return(input_list)

}
