#' @export
#' @import VaRES
#' @import stats
rlse_f <- function(p, alpha, beta) {
    SE <- function(p, alpha, beta) {
        fn = function(x) {
            (VaRES::dF(x,alpha, beta)) * ((VaRES::dF(x,alpha, beta, log = TRUE)) - log(VaRES::pF(p[i],
                alpha, beta)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/VaRES::pF(p, alpha, beta)) * SE(p, alpha, beta)
    ((se_f(alpha, beta) - sinc))/(se_f(alpha, beta))
}

#' @export
rlre_f <- function(p, alpha, beta, delta) {
    RE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dF(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/VaRES::pF(p, alpha, beta))^delta * RE(p,
        alpha, beta, delta))

    (re_f(alpha, beta, delta) - rinc)/re_f(alpha, beta, delta)
}

#' @export
rlhce_f <- function(p, alpha, beta, delta) {
    HC <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dF(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/VaRES::pF(p, alpha, beta))^delta *
        HC(p, alpha, beta, delta)) - 1)
    (hce_f(alpha, beta, delta) - hcinc)/hce_f(alpha, beta, delta)
}

#' @export
rlae_f <- function(p, alpha, beta, delta) {
    AE <- function(p, alpha, beta, delta) {
        fn = function(x) {
            (VaRES::dF(x, alpha, beta))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/VaRES::pF(p, alpha, beta))^delta * AE(p,
        alpha, beta, delta))^(1/delta) - 1)

    (ae_f(alpha, beta, delta) - ainc)/ae_f(alpha, beta, delta)
}
