
#' @export
particle_filter_wrapper <- function(input_list, time1, time2, loss = loss, loss_args){

  output <- particle_filter(input_list$theta, input_list$x, input_list$w, loss = loss, loss_args = loss_args, time1 = time1, time2 = time2)

  size_x <- length(output[[1]]$x)

  output_list <- input_list

  output_list$distance <- sapply(output, function(x){x$distance})

  output_list$x <- matrix(sapply(output, function(x){x$x}), byrow = TRUE, ncol = size_x)

  return(output_list)
}
