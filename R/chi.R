#' @export
#' @import stats

se_chi <- function(n) {
  if (any(n <= 0))
    stop(paste("n must be greater than 0"))

  integrand <- function(x, n) {
    ((stats::dchisq(x, n))*(stats::dchisq(x, n, log = TRUE)))
  }

    fun1=stats::integrate(integrand, lower = 0,upper = Inf, n = n)$value

  return(-1*fun1)
}

#' @export
re_chi <- function(n, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(n <= 0))
    stop(paste("n must be greater than 0"))
  integrand <- function(x, n, delta) {
    (stats::dchisq(x,n))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0, upper = Inf, n = n, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_chi <- function(n, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(n <= 0))
    stop(paste("n must be greater than 0"))
  integrand <- function(x, n, delta) {
    (stats::dchisq(x,n))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0, upper = Inf, n = n, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_chi <- function(n, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(n <= 0))
    stop(paste("n must be greater than 0"))
  integrand <- function(x, n, delta) {
    (stats::dchisq(x,n))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0, upper = Inf, n = n, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


