#' @export
#' @import VaRES
#' @import stats
rlse_gum <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (extraDistr::dgumbel(x,alpha, beta)) * ((extraDistr::dgumbel(x,alpha, beta, log = TRUE)) - log(extraDistr::pgumbel(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower= -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/extraDistr::pgumbel(p, alpha, beta)) * SE(p, alpha, beta)
    ((Se_gum(alpha, beta) - sinc))/(Se_gum(alpha, beta))
}

#' @export
rlre_gum <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dgumbel(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower= -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/extraDistr::pgumbel(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_gum(alpha, beta, delta) - rinc)/re_gum(alpha, beta, delta)
}

#' @export
rlhce_gum <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dgumbel(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower= -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/extraDistr::pgumbel(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_gum(alpha, beta, delta) - hcinc)/hce_gum(alpha, beta, delta)
}

#' @export
rlae_gum <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (extraDistr::dgumbel(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower= -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/extraDistr::pgumbel(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_gum(alpha, beta, delta) - ainc)/ae_gum(alpha, beta, delta)
}
