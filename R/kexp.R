#' @export
#' @import VaRES
#' @import stats

re_kexp <- function(lambda, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
	if (any(b <= 0))
    stop(paste("b must be greater than 0"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))

  integrand <- function(x, lambda, a, b, delta) {
    (VaRES::dkumexp(x, lambda, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda = lambda, a = a, b = b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_kexp <- function(lambda, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(b <= 0))
    stop(paste("b must be greater than 0"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  integrand <- function(x, lambda, a, b, delta) {
    (VaRES::dkumexp(x, lambda, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda = lambda, a = a, b = b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_kexp <- function(lambda, a, b, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(a <= 0))
    stop(paste("a must be greater than 0"))
  if (any(b <= 0))
    stop(paste("b must be greater than 0"))
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  integrand <- function(x, lambda, a, b, delta) {
    (VaRES::dkumexp(x, lambda, a, b))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, lambda = lambda, a = a, b = b, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


