#' @export
#' @import VaRES
#' @import stats
se_ew <- function(a, beta, zeta) {
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))
  integrand <- function(x, a, beta, zeta) {
    ((VaRES::dexpweibull(x, a, beta, zeta))*(VaRES::dexpweibull(x, a, beta, zeta, log = TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = 0,upper = Inf, a = a, beta = beta, zeta = zeta)$value

  return(-1*fun1)
}
#' @export
re_ew <- function(a, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))

  integrand <- function(x, a, beta, zeta, delta) {
    (VaRES::dexpweibull(x, a, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, a = a, beta = beta, zeta = zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_ew <- function(a, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))
  integrand <- function(x, a, beta, zeta, delta) {
    (VaRES::dexpweibull(x, a, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, a = a, beta = beta, zeta = zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_ew <- function(a, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))
  integrand <- function(x, a, beta, zeta, delta) {
    (VaRES::dexpweibull(x, a, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, a = a, beta = beta, zeta = zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


