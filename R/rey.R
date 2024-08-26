#' @export
#' @import stats

se_ray <- function(alpha) {
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  rey<-function(x,alpha){
    2*x*(alpha)*exp(-alpha*x^2)
  }

  rey1<-function(x,alpha){
    log(2)+log(x)+log(alpha)-alpha*x^2
  }

  integrand <- function(x, alpha) {
    ((rey(x, alpha))*(rey1(x, alpha)))
  }

    fun1 <- stats::integrate(integrand, lower = 0,upper = Inf, alpha=alpha)$value

  return(-1*fun1)
}

#' @export
re_ray <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))

  integrand <- function(x, alpha, delta) {
  rey<-function(x,alpha){
2*x*(alpha)*exp(-alpha*x^2)
}
    (rey(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")
  (1/(1 - delta)) * log(fun2(delta))
}

#' @export
hce_ray <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  integrand <- function(x, alpha, delta) {
  rey<-function(x,alpha){
2*x*(alpha)*exp(-alpha*x^2)
}
    (rey(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (1/(2^(1 - delta) - 1)) * (fun2(delta) - 1)
}


#' @export
ae_ray <- function(alpha, delta) {
  if (any(delta == 1))
    stop(paste("delta cannot take exactly 1"))
  if (any(alpha <= 0))
    stop(paste("alpha must be greater than 0"))
  integrand <- function(x, alpha, delta) {
  rey<-function(x,alpha){
2*x*(alpha)*exp(-alpha*x^2)
}
    (rey(x, alpha))^delta
  }
  fun1 <- function(delta) stats::integrate(integrand, lower = 0,
                                           upper = Inf, alpha = alpha, delta = delta)$value
  fun2 <- Vectorize(fun1, "delta")

  (delta/(1 - delta)) * (fun2(delta)^(1/delta) - 1)
}


