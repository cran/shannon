#' @export
#' @import stats
rlse_exp <- function(p, alpha) {
    SE <- function(p, alpha) {
        fn = function(x) {
            (dexp(x, alpha)) * ((dexp(x, alpha, log = TRUE)) - log(pexp(p[i],
                alpha)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/pexp(p, alpha)) * SE(p, alpha)
    ((Se_exp(alpha) - sinc))/(Se_exp(alpha))
}

#' @export
rlre_exp <- function(p, alpha, delta) {
    RE <- function(p, alpha, delta) {
        fn = function(x) {
            (stats::dexp(x, alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/stats::pexp(p, alpha))^delta * RE(p,
        alpha, delta))

    (re_exp(alpha, delta) - rinc)/re_exp(alpha, delta)
}

#' @export
rlhce_exp <- function(p, alpha, delta) {
    HC <- function(p, alpha, delta) {
        fn = function(x) {
            (stats::dexp(x, alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/stats::pexp(p, alpha))^delta *
        HC(p, alpha, delta)) - 1)
    (hce_exp(alpha, delta) - hcinc)/hce_exp(alpha, delta)
}

#' @export
rlae_exp <- function(p, alpha, delta) {
    AE <- function(p, alpha, delta) {
        fn = function(x) {
            (stats::dexp(x, alpha))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/stats::pexp(p, alpha))^delta * AE(p,
        alpha, delta))^(1/delta) - 1)

    (ae_exp(alpha, delta) - ainc)/ae_exp(alpha, delta)
}
