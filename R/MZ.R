#' @title
#' Calculating M-statistics by Stock (1990) and Perron and Ng (1996).
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
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
MZ.statistic <- function(y, l,
                         const = FALSE, trend = FALSE) {
    N <- nrow(y)

    tmp.ADF <- ADF.test(y, const, trend, l, criterion = NULL)

    denom <- 1 - sum(tmp.ADF$beta) + tmp.ADF$alpha
    S.2 <- drop(t(tmp.ADF$residuals) %*% tmp.ADF$residuals) /
        (nrow(tmp.ADF$residuals) - (1 + l)) / denom^2

    sum.y2 <- sum(tmp.ADF$yd[1:(N - 1)]^2)

    mza <- (y[N]^2 / N - S.2) / (2 * sum.y2 / N^2)
    msb <- sqrt(sum.y2 / S.2 / N^2)
    mzt <- mza * msb

    return(
        list(
            mza = mza,
            msb = msb,
            mzt = mzt
        )
    )
}
