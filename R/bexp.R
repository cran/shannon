#' @export
#' @import VaRES
#' @import stats

se_bexp <- function(lambda, alpha, beta) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, lambda, alpha, beta) {
    ((VaRES::dbetaexp(x,lambda, alpha, beta))*(VaRES::dbetaexp(x,lambda, alpha, beta, log=TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = 0,upper = Inf, lambda=lambda, alpha=alpha, beta=beta)$value

  return(-1*fun1)
}

#' @export
re_bexp <- function(lambda, alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))

  integrand <- function(x, lambda, alpha, beta, delta) {
    (VaRES::dbetaexp(x,lambda, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda=lambda, alpha=alpha, beta=beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_bexp <- function(lambda, alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, lambda, alpha, beta, delta) {
    (VaRES::dbetaexp(x,lambda, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda=lambda, alpha=alpha, beta=beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_bexp <- function(lambda, alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, lambda, alpha, beta, delta) {
    (VaRES::dbetaexp(x,lambda, alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda=lambda, alpha=alpha, beta=beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


