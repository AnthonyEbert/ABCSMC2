#' @export
get_parameter <- function(full_list, spread = FALSE){

  TT <- length(full_list)

  Ntheta <- length(full_list[[2]])
  size_theta <- length(full_list[[2]][[1]]$theta)

  theta_array <- array(dim = c(TT, Ntheta, size_theta))
  w_mat <- matrix(nrow = TT, ncol = Ntheta)

  for(tp in 1:TT){
    x_list <- full_list[[tp]]

    theta_array[tp, , ] <- t(sapply(x_list, function(x){x$theta}))

    w_mat[tp, ] <- sapply(x_list, function(x){x$omega}) %>% {. / sum(.)}

  }

  theta_df <- plyr::adply(theta_array, c(1,2,3), .id = c("time", "particle", "parameter")) %>% dplyr::rename(value = V1)

  w_df <- plyr::adply(w_mat, c(1,2))
  names(w_df) <- c("time", "particle", "weight")

  if(!spread){
  theta_df <- dplyr::left_join(theta_df, w_df, by = c("time", "particle")) %>%
    dplyr::group_by(time, parameter, value) %>%
    dplyr::summarise(weight = sum(weight)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time = as.numeric(time))
  } else {
  theta_df <- tidyr::spread(theta_df, parameter, value) %>%
    dplyr::left_join(w_df, by = c("time", "particle")) %>%
    dplyr::group_by(time, `1`, `2`) %>%
    dplyr::summarise(weight = sum(weight)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time = as.numeric(time))
  }

  return(theta_df)
}
