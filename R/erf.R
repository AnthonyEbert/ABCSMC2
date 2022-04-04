
erf <- function(x){2 * pnorm(x * sqrt(2)) - 1}

#' @export
dsn_stan <- function(y, xi, omega, alpha){
  1/(omega * sqrt(2*pi)) * exp(-1/2*((y - xi)/omega)^2) * (1 + erf(alpha*((y - xi)/(omega * sqrt(2)))))
}
