#' @references
#' Perron, Pierre, and Gabriel Rodrı́guez.
#' “GLS Detrending, Efficient Unit Root Tests and Structural Change.”
#' Journal of Econometrics 115, no. 1 (July 1, 2003): 1–27.
#' https://doi.org/10.1016/S0304-4076(03)00090-3.
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
