#' @export
#' @import VaRES
#' @import stats
rlse_kum <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (VaRES::dkum(x,alpha, beta)) * ((VaRES::dkum(x,alpha, beta, log = TRUE)) - log(VaRES::pkum(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/VaRES::pkum(p, alpha, beta)) * SE(p, alpha, beta)
    ((se_kum(alpha, beta) - sinc))/(se_kum(alpha, beta))
}

#' @export
rlre_kum <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dkum(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/VaRES::pkum(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_kum(alpha, beta, delta) - rinc)/re_kum(alpha, beta, delta)
}

#' @export
rlhce_kum <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dkum(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/VaRES::pkum(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_kum(alpha, beta, delta) - hcinc)/hce_kum(alpha, beta, delta)
}

#' @export
rlae_kum <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dkum(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/VaRES::pkum(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_kum(alpha, beta, delta) - ainc)/ae_kum(alpha, beta, delta)
}
