#' @export
#' @import VaRES
#' @import stats

se_kumnorm <- function(mu, sigma, a, b) {
  if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(b <= 0))
    stop(paste("b must be greater than 0"))

  integrand <- function(x, mu, sigma, a, b) {
    (VaRES::dkumnormal(x, mu, sigma, a, b))^0.99999
  }
  fun1 <- stats::integrate(integrand, lower = -Inf,
                          upper = Inf, mu=mu, sigma=sigma, a=a, b=b)$value
  z<-(1/(1 - 0.99999)) * log(fun1)
  return(round(z,4))
}

#' @export
re_kumnorm <- function(mu, sigma, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
	if (any(b <= 0))
    stop(paste("b must be greater than 0"))

  integrand <- function(x, mu, sigma, a, b, delta) {
    (VaRES::dkumnormal(x, mu, sigma, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, mu=mu, sigma=sigma, a=a, b=b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_kumnorm <- function(mu, sigma, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
	if (any(b <= 0))
    stop(paste("b must be greater than 0"))
  integrand <- function(x, mu, sigma, a, b, delta) {
    (VaRES::dkumnormal(x, mu, sigma, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, mu=mu, sigma=sigma, a=a, b=b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_kumnorm <- function(mu, sigma, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
	if (any(b <= 0))
    stop(paste("b must be greater than 0"))
  integrand <- function(x, mu, sigma, a, b, delta) {
    (VaRES::dkumnormal(x, mu, sigma, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, mu=mu, sigma=sigma, a=a, b=b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


