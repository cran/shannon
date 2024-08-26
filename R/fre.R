#' @export
#' @import stats


se_fre <- function(alpha, beta, zeta) {
 if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))

  integrand <- function(x, alpha, beta, zeta) {
    ((extraDistr::dfrechet(x, alpha, beta, zeta))*(extraDistr::dfrechet(x, alpha, beta, zeta, log = TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = beta, upper = Inf, alpha = alpha, beta=beta, zeta=zeta)$value

  return(-1*fun1)
}

#' @export
re_fre <- function(alpha, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))

  integrand <- function(x, alpha, beta, zeta, delta) {
    (extraDistr::dfrechet(x, alpha, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = beta, upper = Inf, alpha = alpha, beta=beta, zeta=zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_fre <- function(alpha, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))

  integrand <- function(x, alpha, beta, zeta, delta) {
    (extraDistr::dfrechet(x, alpha, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = beta, upper = Inf, alpha = alpha, beta=beta, zeta=zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_fre <- function(alpha, beta, zeta, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  if (any(zeta <= 0))
    stop(paste("zeta must be greater than 0"))

  integrand <- function(x, alpha, beta, zeta, delta) {
    (extraDistr::dfrechet(x, alpha, beta, zeta))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = beta, upper = Inf, alpha = alpha, beta=beta, zeta=zeta, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


