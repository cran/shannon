#' @export
#' @import stats
rlse_norm <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (stats::dnorm(x,alpha, beta)) * ((stats::dnorm(x,alpha, beta, log = TRUE)) - log(stats::pnorm(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/stats::pnorm(p, alpha, beta)) * SE(p, alpha, beta)
    ((se_norm(alpha, beta) - sinc))/(se_norm(alpha, beta))
}

#' @export
rlre_norm <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dnorm(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/stats::pnorm(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_norm(alpha, beta, delta) - rinc)/re_norm(alpha, beta, delta)
}

#' @export
rlhce_norm <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dnorm(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/stats::pnorm(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_norm(alpha, beta, delta) - hcinc)/hce_norm(alpha, beta, delta)
}

#' @export
rlae_norm <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (stats::dnorm(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/stats::pnorm(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_norm(alpha, beta, delta) - ainc)/ae_norm(alpha, beta, delta)
}
