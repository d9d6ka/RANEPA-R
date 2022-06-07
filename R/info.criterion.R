#' @title
#' Information criterions
#'
#' @details
#' Calculating the value of the following informational criterions:
#' \itemize{
#' \item Akaike,
#' \item Schwarz (Bayesian),
#' \item Hannan-Quinn,
#' \item Liu et al.}
#'
#' @param resid Input residuals.
#' @param extra Number of coefficients etc.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param modification Whether the unit-root test modificaton is needed.
#' @param alpha The coefficient \eqn{\alpha} of \eqn{y_{t-1}} in ADF model.
#' @param y The vector of \eqn{y_{t-1}} in ADF model.
#'
#' @return
#' The value of the selected informational criterion.
#'
#' @export
info.criterion <- function(resid, extra,
                           modification = FALSE,
                           alpha = 0, y = NULL) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)

    N <- nrow(resid)

    if (modification) {
        s2 <- drop(t(resid) %*% resid) / (N - 1)
        tau <- (alpha ^ 2) * drop(t(y) %*% y) / s2
    } else {
        tau <- 0
    }

    log.RSS <- log(drop(t(resid) %*% resid) / N)

    return(
        list(
            aic = log.RSS + 2 * (tau + extra) / N,
            bic = log.RSS + (tau + extra) * log(N) / N,
            hq = log.RSS + 2 * (tau + extra) * log(log(N)) / N,
            lwz = log.RSS + 0.299 * (tau + extra) * (log(N))^2.1
        )
    )
}

#' @importFrom zeallot %<-%
lag.selection <- function(y, x, max.lag, criterion = "aic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)
    k <- ncol(x)

    tmp.y <- y[(1 + max.lag):N, , drop = FALSE]
    tmp.x <- x[(1 + max.lag):N, , drop = FALSE]

    for (l in 1:max.lag) {
        tmp.x <- cbind(
            tmp.x,
            lagn(y, l)[(1 + max.lag):N, , drop = FALSE]
        )
    }

    c(., resid, ., .) %<-% OLS(tmp.y, tmp.x[, 1:k, drop = FALSE])

    res.IC <- log(drop(t(resid) %*% resid) / (N - max.lag))
    res.lag <- 0

    for (l in 1:max.lag) {
        c(., resid, ., .) %<-% OLS(tmp.y, tmp.x[, 1:(k + l), drop = FALSE])
        temp.IC <- info.criterion(resid, l)[[criterion]]

        if (temp.IC < res.IC) {
            res.IC <- temp.IC
            res.lag <- l
        }
    }

    return(res.lag)
}

#' @importFrom zeallot %<-%
ADF.lag.selection <- function(y, const = FALSE, trend = FALSE, max.lag,
                              criterion = "aic", modification = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)
    pos <- ifelse(const && trend, 3, ifelse(const, 2, 1))

    d.y <- c(NA, diff(y))
    x <- cbind(
        if (const) rep(1, N) else NULL,
        if (trend) 1:N else NULL,
        lagn(y, 1)
    )

    for (l in 1:max.lag) {
        x <- cbind(x, lagn(d.y, l))
    }

    tmp.y <- d.y[(2 + max.lag):N, , drop = FALSE]
    tmp.x <- x[(2 + max.lag):N, , drop = FALSE]

    c(., resid, ., .) %<-% OLS(tmp.y, tmp.x[, 1:pos, drop = FALSE])

    res.IC <- log(drop(t(resid) %*% resid) / (N - max.lag))
    res.lag <- 0

    for (l in 1:max.lag) {
        c(beta, resid, ., .) %<-% OLS(tmp.y, tmp.x[, 1:(pos + l), drop = FALSE])
        temp.IC <- info.criterion(
            resid, l,
            modification, beta[pos], tmp.x[, pos, drop = FALSE]
        )[[criterion]]

        if (temp.IC < res.IC) {
            res.IC <- temp.IC
            res.lag <- l
        }
    }

    return(res.lag)
}
