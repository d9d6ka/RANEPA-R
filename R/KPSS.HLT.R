#' Unit root testing procedure under a single structural break.
#'
#' @param y The series of interest.
#' @param const Whether a constant should be included.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and A. M. Robert Taylor.
#' “Unit Root Testing under a Local Break in Trend.”
#' Journal of Econometrics 167, no. 1 (2012): 140–67.
#'
#' @importFrom zeallot %<-%
#' @export
KPSS.HLT <- function(y, const = FALSE, trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    if (!const) {
        m.ksi <- 0.853
    } else {
        m.ksi <- 1.052
    }

    bartlett.lag <- trunc(4 * ((N / 100)^(1 / 4)))

    d.y <- diff(y)

    first.lag <- trunc(trim * N)
    last.lag <- trunc((1 - trim) * N)

    t0 <- -Inf
    t1 <- -Inf

    tb0 <- NA
    tb1 <- NA

    for (tb in first.lag:last.lag) {
        tau <- tb / N

        DU <- c(rep(0, tb), rep(1, N - tb))
        DT <- DU * (1:N - tb)

        x <- cbind(
            rep(1, N),
            1:N,
            if (const) {
                DU
            } else {
                NULL
            },
            DT
        )

        c(beta, resid, ., .) %<-% OLS(y, x)

        lr.var.y <- lr.var.bartlett(resid, bartlett.lag)
        inv.xx <- qr.solve(t(x) %*% x)

        t0.stat <-
            abs(beta[ncol(x)] / sqrt(lr.var.y * inv.xx[ncol(x), ncol(x)]))

        x <- cbind(
            rep(1, N - 1),
            if (const) diff(DU) else NULL,
            DU[2:N]
        )

        c(beta, resid, ., .) %<-% OLS(d.y, x)

        lr.var.dy <- lr.var.bartlett(resid, bartlett.lag)
        inv.xx <- qr.solve(t(x) %*% x)

        t1.stat <-
            abs(beta[ncol(x)] / sqrt(lr.var.dy * inv.xx[ncol(x), ncol(x)]))

        if (t0.stat > t0) {
            t0 <- t0.stat
            tb0 <- tb
        }
        if (t1.stat > t1) {
            t1 <- t1.stat
            tb1 <- tb
        }
    }

    DU0 <- c(rep(0, tb0), rep(1, N - tb0))
    DT0 <- DU0 * (1:N - tb0)

    x <- cbind(
        rep(1, N),
        1:N,
        if (const) DU0 else NULL,
        DT0
    )

    c(beta, resid, ., .) %<-% OLS(y, x)

    lr.var.y <- lr.var.bartlett(resid, bartlett.lag)
    kpss.y <- KPSS(resid, lr.var.y)

    DU1 <- c(rep(0, tb1), rep(1, N - tb1))

    x <- cbind(
        rep(1, N - 1),
        if (const) diff(DU1) else NULL,
        DU1[2:N]
    )

    c(beta, resid, ., .) %<-% OLS(d.y, x)

    lr.var.dy <- lr.var.bartlett(resid, bartlett.lag)
    kpss.dy <- KPSS(resid, lr.var.dy)

    lambda.kpss <- exp(-((500 * kpss.y * kpss.dy)^2))
    t.lambda.kpss <- lambda.kpss * t0 + m.ksi * (1 - lambda.kpss) * t1

    return(t.lambda.kpss)
}
