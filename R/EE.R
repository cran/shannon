#' @export
#' @import VaRES
#' @import stats
se_ee <- function(alpha, beta) {
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta) {
    ((VaRES::dexpexp(x, alpha, beta))*(VaRES::dexpexp(x, alpha, beta, log = TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = 0,upper = Inf, alpha=alpha, beta=beta)$value

  return(-1*fun1)
}

#' @export
re_ee <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (VaRES::dexpexp(x, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_ee <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (VaRES::dexpexp(x, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_ee <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta, delta) {
    (VaRES::dexpexp(x, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


