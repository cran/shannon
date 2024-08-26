#' @export
#' @import VaRES
#' @import stats

Se_gum <- function(alpha, beta) {
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))

  integrand <- function(x, alpha, beta) {
    (extraDistr::dgumbel(x,alpha, beta))^0.99999
  }

  fun1 <- stats::integrate(integrand, lower=-Inf,
                           upper = Inf, alpha = alpha, beta = beta)$value

  z<-(1/(1 - 0.99999)) * log(fun1)
  return(round(z,4))
}


#' @export
re_gum <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))

  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dgumbel(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_gum <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dgumbel(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_gum <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dgumbel(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


