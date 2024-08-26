#' @export
#' @import stats

se_logis <- function(mu, sigma) {
  if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  integrand <- function(x, mu, sigma) {
    ((stats::dlogis(x, mu, sigma))*(stats::dlogis(x, mu, sigma, log = TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = -Inf,upper = Inf, mu=mu, sigma=sigma)$value

  return(-1*fun1)
}

#' @export
re_logis <- function(mu, sigma, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))

  integrand <- function(x, mu, sigma, delta) {
    (stats::dlogis(x,mu, sigma))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,upper = Inf, mu = mu, sigma = sigma, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_logis <- function(mu, sigma, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  integrand <- function(x, mu, sigma, delta) {
    (stats::dlogis(x,mu, sigma))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,upper = Inf, mu = mu, sigma = sigma, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_logis <- function(mu, sigma, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0"))
  integrand <- function(x, mu, sigma, delta) {
    (stats::dlogis(x,mu, sigma))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,upper = Inf, mu = mu, sigma = sigma, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


