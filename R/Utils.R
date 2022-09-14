#' @title
#' Produce a vector lagged backward of forward
#'
#' @param x Initial vector.
#' @param i Size of lag (lead if negative).
#' @param na Value to fill missing observations, `NA` by default.
#'
#' @return Lagged or leaded vector.
#'
#' @keywords internal
lagn <- function(x,
                 i,
                 na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    n.obs <- nrow(x)
    k <- ncol(x)
    if (i > 0) {
        return(
            rbind(
                matrix(data = na, nrow = i, ncol = k),
                x[1:(n.obs - i), , drop = FALSE]
            )
        )
    } else {
        return(
            rbind(
                x[(1 + abs(i)):n.obs, , drop = FALSE],
                matrix(data = na, nrow = abs(i), ncol = k)
            )
        )
    }
}


#' @title
#' Calculating ACF values
#'
#' @description
#' A simple auxiliary function providing the estimates of autocorrelations of
#' orders from \eqn{0} to \eqn{N - 1}.
#'
#' @param y An input time series of interest.
#'
#' @return The vector of ACF values of orders from 0 to \eqn{N - 1}.
#' Due to the R's way of indexing the 1-st element is the autocorrelation of
#' order 0, and the N-th value is the autcorrelation of order \eqn{N - 1}.
#'
#' @keywords internal
ACF <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    mean.y <- mean(y)

    result <- rep(0, n.obs)

    for (i in 0:(n.obs - 1)) {
        result[i + 1] <-
            t(y[1:(n.obs - i)] - mean.y) %*% (y[(1 + i):n.obs] - mean.y) / n.obs
    }

    return(result)
}


#' @title
#' Auxiliary function returning KPSS statistic value.
#'
#' @param resids The series of residuals.
#' @param variance The value of the long-run variance.
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
#' @param y A series of interest.
#' @param l Number of lags for inner ADF test.
#' @param const,trend Whether a constand and trend are to be included.
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
