#' @export
#' @import stats
Se_exp <- function(alpha) {
  integrand <- function(x, alpha) {
    ((stats::dexp(x, alpha))*(stats::dexp(x, alpha, log=TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = 0,upper = Inf, alpha=alpha)$value

  return(-1*fun1)
}

#' @export
re_exp <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))

  integrand <- function(x, alpha, delta) {
    (stats::dexp(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_exp <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  integrand <- function(x, alpha, delta) {
    (stats::dexp(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_exp <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  integrand <- function(x, alpha, delta) {
    (stats::dexp(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


