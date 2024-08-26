#' @export
#' @import VaRES
#' @import stats
rlse_ig <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (extraDistr::dinvgamma(x,alpha, beta)) * ((extraDistr::dinvgamma(x,alpha, beta, log = TRUE)) - log(extraDistr::pinvgamma(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/extraDistr::pinvgamma(p, alpha, beta)) * SE(p, alpha, beta)
    ((se_ig(alpha, beta) - sinc))/(se_ig(alpha, beta))
}

#' @export
rlre_ig <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dinvgamma(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/extraDistr::pinvgamma(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_ig(alpha, beta, delta) - rinc)/re_ig(alpha, beta, delta)
}

#' @export
rlhce_ig <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dinvgamma(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/extraDistr::pinvgamma(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_ig(alpha, beta, delta) - hcinc)/hce_ig(alpha, beta, delta)
}

#' @export
rlae_ig <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dinvgamma(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/extraDistr::pinvgamma(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_ig(alpha, beta, delta) - ainc)/ae_ig(alpha, beta, delta)
}
