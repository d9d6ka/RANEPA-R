#' @importFrom zeallot %<-%
#'
#' @export
robust.tests <- function(y,
                         const = FALSE, trend = FALSE,
                         trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    cv.DF.GLS.m <- -1.93078
    cv.DF.OLS.m <- -2.85706

    cv.DF.GLS.t <- -2.84317
	cv.MDF.GLS <- -3.84632
	cv.DF.OLS.t <- -3.39735
	cv.MDF.OLS <- -4.23
	cv.MDF.t <- -4.25

    if (const && trend) {
        cv.MDF.OLS <- -4.21172
        cv.MDF.t <- -4.34468
    }

    if (const && !trend) cv.HLT <- 2.563
    if (const && trend) cv.HLT <- 3.162

	sap.ur <- 1.1009;
	sap.ur1 <- 1.0514;
	sap.ur2 <- 1.0317;
	sap.ur3 <- 1.0364;

	sap.cv.ur.HLT.k0 <- 1.0154
	sap.cv.A.HLT.k0 <- 1.00

	sap.cv.ur.PY.k0 <- 1.0159
	sap.cv.A.PY.k0 <- 1.00

    N <- nrow(y)

    max.lag <- trunc(12 * (N / 100)^(1 / 4))

    first.break <- trunc(trim * N)
    last.break <- trunc((1 - trim) * N)

    tb <- segs.GLS.1.break(
        y, const, trend,
        first.break, last.break,
        trim
    )
    tau <- tb / N

    DU <- c(rep(0, tb), rep(1, N - tb))
    DT <- DU * (1:N - tb)

    x <- cbind(
        rep(1, N),
        if (const) DU else NULL,
        1:N,
        if (trend) DT else NULL
    )

    ## Mean case
    c(beta.OLS.m, resid.OLS.m, ., .) %<-% OLS(y, x[, 1])
    k.m <- ADF.lag.selection(
        y = y, max.lag = max.lag,
        criterion = "aic", modification = TRUE
    )

    res.ADF <- ADF.test(y, const = FALSE, trend = FALSE,
                        max.lag = k.m)

    resid.GLS.m <- GLS(y, x[, 1], -7)

}
