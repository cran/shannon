#' @export
#' @import stats
rlse_wei <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (stats::dweibull(x,alpha, beta)) * ((stats::dweibull(x,alpha, beta, log = TRUE)) - log(stats::pweibull(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/stats::pweibull(p, alpha, beta)) * SE(p, alpha, beta)
    ((se_wei(alpha, beta) - sinc))/(se_wei(alpha, beta))
}

#' @export
rlre_wei <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dweibull(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/stats::pweibull(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_wei(alpha, beta, delta) - rinc)/re_wei(alpha, beta, delta)
}

#' @export
rlhce_wei <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dweibull(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/stats::pweibull(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_wei(alpha, beta, delta) - hcinc)/hce_wei(alpha, beta, delta)
}

#' @export
rlae_wei <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dweibull(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/stats::pweibull(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_wei(alpha, beta, delta) - ainc)/ae_wei(alpha, beta, delta)
}
