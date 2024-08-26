#' @export
#' @import stats
se_st <- function(v) {
  if (any(v <= 0))
    stop(paste("v must be greater than 0"))
  integrand <- function(x, v) {
    ((stats::dt(x, v))*(stats::dt(x, v, log=TRUE)))
  }

    fun1 <- stats::integrate(integrand, lower = -Inf,upper = Inf, v=v)$value

  return(-1*fun1)
}

#' @export
re_st <- function(v, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(v <= 0))
    stop(paste("v must be greater than 0"))

  integrand <- function(x, v, delta) {
    (stats::dt(x, v))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, v = v, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_st <- function(v, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(v <= 0))
    stop(paste("v must be greater than 0"))
  integrand <- function(x, v, delta) {
    (stats::dt(x, v))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, v = v, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}

#' @export
ae_st <- function(v, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(v <= 0))
    stop(paste("v must be greater than 0"))
  integrand <- function(x, v, delta) {
    (stats::dt(x, v))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = -Inf,
                                           upper = Inf, v = v, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


