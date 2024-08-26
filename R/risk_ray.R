#' @export
#' @import stats

rlse_ray <- function(p, alpha) {
    drey<-function(x,alpha) {
    2*x*(alpha)*exp(-alpha*x^2)
  }
  log_rey<-function(x,alpha){
    log(2)+log(x)+log(alpha)-alpha*x^2
  }
  prey<-function(p, alpha) {
    1-exp(-1*alpha*p^2)
}

    SE <- function(p, alpha) {
        fn = function(x) {
          (drey(x,alpha)) * ((log_rey(x, alpha)) - (log(1-exp(-alpha*p[i]^2))))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- (-1 * (1/prey(p, alpha)) * SE(p, alpha))
    ((se_ray(alpha) - sinc))/(se_ray(alpha))
}

#' @export
rlre_ray <- function(p, alpha, delta) {

 drey<-function(x,alpha) {
    2*x*(alpha)*exp(-alpha*x^2)
  }
  prey<-function(p, alpha) {
    1-exp(-alpha*p^2)
  }

    RE <- function(p, alpha, delta) {
        fn = function(x) {
            (drey(x,alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/prey(p, alpha))^delta * RE(p,
        alpha, delta))

    (re_ray(alpha, delta) - rinc)/re_ray(alpha, delta)
}

#' @export
rlhce_ray <- function(p, alpha, delta) {
drey<-function(x,alpha) {
    2*x*(alpha)*exp(-alpha*x^2)
  }
  prey<-function(p, alpha) {
    1-exp(-alpha*p^2)
  }
    HC <- function(p, alpha, delta) {
        fn = function(x) {
            (drey(x,alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/prey(p, alpha))^delta *
        HC(p, alpha, delta)) - 1)
    (hce_ray(alpha, delta) - hcinc)/hce_ray(alpha, delta)
}

#' @export
rlae_ray <- function(p, alpha, delta) {
drey<-function(x,alpha) {
    2*x*(alpha)*exp(-alpha*x^2)
  }
  prey<-function(p, alpha) {
    1-exp(-alpha*p^2)
  }
    AE <- function(p, alpha, delta) {
        fn = function(x) {
            (drey(x,alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/prey(p, alpha))^delta * AE(p,
        alpha, delta))^(1/delta) - 1)

    (ae_ray(alpha, delta) - ainc)/ae_ray(alpha, delta)
}
