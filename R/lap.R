#' @export
#' @import VaRES
#' @import stats

Se_lap <- function(alpha, beta) {
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta) {
    ((extraDistr::dlaplace(x, alpha, beta))*(extraDistr::dlaplace(x, alpha, beta, log=TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower=-Inf,upper=Inf, alpha=alpha, beta=beta)$value

  return(-1*fun1)
}

#' @export
re_lap <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))

  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dlaplace(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_lap <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dlaplace(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_lap <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (extraDistr::dlaplace(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=-Inf,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


