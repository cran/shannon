#' @export
#' @import VaRES
#' @import stats

se_burr <- function(k, c) {

  if (any(k <= 0))
    stop(paste("k must be greater than 0"))
  if (any(c <= 0))
    stop(paste("c must be greater than 0"))
  integrand <- function(x, k, c) {
    ((VaRES::dburr7(x, k, c))*(VaRES::dburr7(x, k, c, log = TRUE)))
  }

    fun1=stats::integrate(integrand, lower=0,upper=Inf, k=k, c=c)$value

  return(-1*fun1)
}

#' @export
re_burr <- function(k, c, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(k <= 0))
    stop(paste("k must be greater than 0"))
  if (any(c <= 0))
    stop(paste("c must be greater than 0"))

  integrand <- function(x, k, c, delta) {
    (VaRES::dburr7(x, k, c))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=0,
                                           upper = Inf, k = k, c = c, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_burr <- function(k, c, delta) {
 if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(k <= 0))
    stop(paste("k must be greater than 0"))
  if (any(c <= 0))
    stop(paste("c must be greater than 0"))
  integrand <- function(x, k, c, delta) {
    (VaRES::dburr7(x, k, c))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=0,
                                           upper = Inf, k = k, c = c, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_burr <- function(k, c, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(k <= 0))
    stop(paste("k must be greater than 0"))
  if (any(c <= 0))
    stop(paste("c must be greater than 0"))
  integrand <- function(x, k, c, delta) {
    (VaRES::dburr7(x, k, c))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower=0,
                                           upper = Inf, k = k, c = c, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


