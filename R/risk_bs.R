#' @export
#' @import stats
rlse_bs <- function(p, v) {
    SE <- function(p, v) {
        fn = function(x) {
            (VaRES::dBS(x, v)) * ((VaRES::dBS(x, v, log = TRUE)) - log(VaRES::pBS(p[i],
                v)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/VaRES::pBS(p, v)) * SE(p, v)
    ((se_bs(v) - sinc))/(se_bs(v))
}

#' @export
rlre_bs <- function(p, v, delta) {
    RE <- function(p, v, delta) {
        fn = function(x) {
            (VaRES::dBS(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/VaRES::pBS(p, v))^delta * RE(p,
        v, delta))

    (re_bs(v, delta) - rinc)/re_bs(v, delta)
}

#' @export
rlhce_bs <- function(p, v, delta) {
    HC <- function(p, v, delta) {
        fn = function(x) {
            (VaRES::dBS(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/VaRES::pBS(p, v))^delta *
        HC(p, v, delta)) - 1)
    (hce_bs(v, delta) - hcinc)/hce_bs(v, delta)
}

#' @export
rlae_bs <- function(p, v, delta) {
    AE <- function(p, v, delta) {
        fn = function(x) {
            (VaRES::dBS(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = 0, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/VaRES::pBS(p, v))^delta * AE(p,
        v, delta))^(1/delta) - 1)

    (ae_bs(v, delta) - ainc)/ae_bs(v, delta)
}
