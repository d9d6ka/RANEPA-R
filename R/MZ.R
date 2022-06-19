MZ.statistic <- function(y, l,
                         const = FALSE, trend = FALSE) {
    N <- nrow(y)

    tmp.ADF <- ADF.test(y, const, trend, l, criterion = NULL)

    denom <- 1 - sum(tmp.ADF$beta) + tmp.ADF$alpha
    S.2 <- drop(t(tmp.ADF$residuals) %*% tmp.ADF$residuals) /
        nrow(tmp.ADF$residuals) / denom

    sum.y2 <- sum(tmp.ADF$yd[1:(N - 1)]^2)

    mza <- (y[N] / N - S.2) / (2 * sum.y2 / N^2)
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
