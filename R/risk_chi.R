#' @export
#' @import stats
rlse_chi <- function(p, n) {
    SE <- function(p, n) {
        fn = function(x) {
            (stats::dchisq(x, n)) * ((stats::dchisq(x, n, log = TRUE)) - log(stats::pchisq(p[i],
                n)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/stats::pchisq(p, n)) * SE(p, n)
    ((se_chi(n) - sinc))/(se_chi(n))
}

#' @export
rlre_chi <- function(p, n, delta) {
    RE <- function(p, n, delta) {
        fn = function(x) {
            (stats::dchisq(x, n))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/stats::pchisq(p, n))^delta * RE(p,
        n, delta))

    (re_chi(n, delta) - rinc)/re_chi(n, delta)
}

#' @export
rlhce_chi <- function(p, n, delta) {
    HC <- function(p, n, delta) {
        fn = function(x) {
            (stats::dchisq(x, n))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/stats::pchisq(p, n))^delta *
        HC(p, n, delta)) - 1)
    (hce_chi(n, delta) - hcinc)/hce_chi(n, delta)
}

#' @export
rlae_chi <- function(p, n, delta) {
    AE <- function(p, n, delta) {
        fn = function(x) {
            (stats::dchisq(x, n))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/stats::pchisq(p, n))^delta * AE(p,
        n, delta))^(1/delta) - 1)

    (ae_chi(n, delta) - ainc)/ae_chi(n, delta)
}
