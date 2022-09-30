#' @title
#' Produce a vector lagged backward of forward
#'
#' @param x An initial vector.
#' @param i Size of lag (lead if negative).
#' @param na A value to fill missing observations, `NA` by default.
#'
#' @return Lagged or leaded vector.
#'
#' @keywords internal
lagn <- function(x,
                 i,
                 na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    n.var <- ncol(x)
    if (i > 0) {
        return(
            rbind(
                matrix(data = na, nrow = i, ncol = n.var),
                x[1:(n.obs - i), , drop = FALSE]
            )
        )
    } else {
        return(
            rbind(
                x[(1 + abs(i)):n.obs, , drop = FALSE],
                matrix(data = na, nrow = abs(i), ncol = n.var)
            )
        )
    }
}


#' @title
#' Produce a vector or matrix of differences, keeping initial length
#'
#' @param x An initial vector.
#' @param lag Size of lag.
#' @param difference Order of differentiating.
#' @param na A value to fill missing observations, `NA` by default.
#'
#' @return Vector or matrix of differences.
#'
#' @keywords internal
diffn <- function(x,
                  lag = 1,
                  differences = 1,
                  na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    n.var <- ncol(x)
    tmp.diff <- diff(x, lag = lag, differences = differences)
    return(
        rbind(
            matrix(
                data = na,
                nrow = n.obs - nrow(tmp.diff),
                ncol = n.var
            ),
            tmp.diff
        )
    )
}


#' @title
#' Auxiliary function returning KPSS statistic value.
#'
#' @param resids A series of residuals.
#' @param variance A value of the long-run variance.
#'
#' @keywords internal
KPSS <- function(resids,
                 variance) {
    if (!is.matrix(resids)) resids <- as.matrix(resids)

    n.obs <- nrow(resids)

    S.t <- apply(resids, 2, cumsum)

    return(drop(t(S.t) %*% S.t) / (n.obs^2 * variance))
}


#' @title
#' Calculating M-statistics by Stock (1990) and Perron and Ng (1996).
#'
#' @param y A time series of interest.
#' @param l Number of lags for inner ADF test.
#' @param const,trend Whether a constant and trend are to be included.
#'
#' @return List of values of \eqn{MZ_\alpha}, \eqn{MZ_t} and \eqn{MSB}
#' statistics.
#'
#' @references
#' Perron, Pierre, and Serena Ng.
#' “Useful Modifications to Some Unit Root Tests with Dependent Errors
#' and Their Local Asymptotic Properties.”
#' The Review of Economic Studies 63, no. 3 (July 1, 1996): 435–63.
#' https://doi.org/10.2307/2297890.
#'
#' Stock, James H.
#' “A Class of Tests for Integration and Cointegration.”
#' Kennedy School of Government, Harvard University, 1990.
#'
#' @keywords internal
MZ.statistic <- function(y,
                         l,
                         const = FALSE,
                         trend = FALSE) {
    n.obs <- nrow(y)

    tmp.ADF <- ADF.test(y, const, trend, l, criterion = NULL)

    denom <- 1 - sum(tmp.ADF$beta) + tmp.ADF$alpha
    S.2 <- drop(t(tmp.ADF$residuals) %*% tmp.ADF$residuals) /
        (nrow(tmp.ADF$residuals) - (1 + l)) / denom^2

    sum.y2 <- sum(tmp.ADF$yd[1:(n.obs - 1)]^2)

    mza <- (y[n.obs]^2 / n.obs - S.2) / (2 * sum.y2 / n.obs^2)
    msb <- sqrt(sum.y2 / S.2 / n.obs^2)
    mzt <- mza * msb

    return(
        list(
            mza = mza,
            msb = msb,
            mzt = mzt
        )
    )
}
