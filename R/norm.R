#' @export
#' @import stats

se_norm <- function(alpha, beta) {
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
  integrand <- function(x, alpha, beta) {
    ((stats::dnorm(x,alpha, beta))*(stats::dnorm(x, alpha, beta, log=TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = -Inf, upper = Inf, alpha=alpha, beta=beta)$value

  return(-1*fun1)
}

#' @export
re_norm <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))

  integrand <- function(x, alpha, beta, delta) {
    (stats::dnorm(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf, upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_norm <- function(alpha, beta, delta) {
 if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
 integrand <- function(x, alpha, beta, delta) {
    (stats::dnorm(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf, upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_norm <- function(alpha, beta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
	if (any(beta <= 0))
    stop(paste("beta must be greater than 0"))
 integrand <- function(x, alpha, beta, delta) {
    (stats::dnorm(x,alpha, beta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf, upper = Inf, alpha = alpha, beta = beta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


