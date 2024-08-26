#' @export
#' @import stats
rlse_st <- function(p, v) {
    SE <- function(p, v) {
        fn = function(x) {
            (stats::dt(x, v)) * ((stats::dt(x, v, log = TRUE)) - log(stats::pt(p[i],
                v)))
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    sinc <- -1 * (1/stats::pt(p, v)) * SE(p, v)
    ((se_st(v) - sinc))/(se_st(v))
}

#' @export
rlre_st <- function(p, v, delta) {
    RE <- function(p, v, delta) {
        fn = function(x) {
            (stats::dt(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }

    rinc <- (1/(1 - delta)) * log((1/stats::pt(p, v))^delta * RE(p,
        v, delta))

    (re_st(v, delta) - rinc)/re_st(v, delta)
}

#' @export
rlhce_st <- function(p, v, delta) {
    HC <- function(p, v, delta) {
        fn = function(x) {
            (stats::dt(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    hcinc <- (1/(2^(1 - delta) - 1)) * (((1/stats::pt(p, v))^delta *
        HC(p, v, delta)) - 1)
    (hce_st(v, delta) - hcinc)/hce_st(v, delta)
}

#' @export
rlae_st <- function(p, v, delta) {
    AE <- function(p, v, delta) {
        fn = function(x) {
            (stats::dt(x, v))^delta
        }
        ES = p
        for (i in 1:length(p)) {
            ES[i] = (stats::integrate(fn, lower = -Inf, upper = p[i], stop.on.error = FALSE)$value)
        }
        return(ES)
    }
    ainc <- (delta/(1 - delta)) * (((1/stats::pt(p, v))^delta * AE(p,
        v, delta))^(1/delta) - 1)

    (ae_st(v, delta) - ainc)/ae_st(v, delta)
}
