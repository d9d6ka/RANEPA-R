#' @title
#' Unit root testing procedure under a single structural break.
#'
#' @param y A series of interest.
#' @param const Whether a constant should be included.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#'
#' @return The value of test statistic.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and A. M. Robert Taylor.
#' “Unit Root Testing under a Local Break in Trend.”
#' Journal of Econometrics 167, no. 1 (2012): 140–67.
#'
#' @export
KPSS.HLT <- function(y, const = FALSE, trim = 0.15) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    if (!const) {
        m.ksi <- 0.853
    } else {
        m.ksi <- 1.052
    }

    d.y <- diff(y)

    first.lag <- trunc(trim * n.obs)
    last.lag <- trunc((1 - trim) * n.obs)

    t0 <- -Inf
    t1 <- -Inf

    tb0 <- NA
    tb1 <- NA

    for (tb in first.lag:last.lag) {
        DU <- c(rep(0, tb), rep(1, n.obs - tb))
        DT <- DU * (1:n.obs - tb)

        x <- cbind(
            rep(1, n.obs),
            1:n.obs,
            if (const) {
                DU
            } else {
                NULL
            },
            DT
        )

        tmp.OLS <- OLS(y, x)

        lr.var.y <- lr.var.bartlett(tmp.OLS$residuals)
        inv.xx <- qr.solve(t(x) %*% x)

        t0.stat <- abs(tmp.OLS$beta[ncol(x)] /
                       sqrt(lr.var.y * inv.xx[ncol(x), ncol(x)]))

        x <- cbind(
            rep(1, n.obs - 1),
            if (const) diff(DU) else NULL,
            DU[2:n.obs]
        )

        tmp.OLS <- OLS(d.y, x)

        lr.var.dy <- lr.var.bartlett(tmp.OLS$residuals)
        inv.xx <- qr.solve(t(x) %*% x)

        t1.stat <- abs(tmp.OLS$beta[ncol(x)] /
                       sqrt(lr.var.dy * inv.xx[ncol(x), ncol(x)]))

        if (t0.stat > t0) {
            t0 <- t0.stat
            tb0 <- tb
        }
        if (t1.stat > t1) {
            t1 <- t1.stat
            tb1 <- tb
        }
    }

    DU0 <- c(rep(0, tb0), rep(1, n.obs - tb0))
    DT0 <- DU0 * (1:n.obs - tb0)

    x <- cbind(
        rep(1, n.obs),
        1:n.obs,
        if (const) DU0 else NULL,
        DT0
    )

    tmp.OLS <- OLS(y, x)

    lr.var.y <- lr.var.bartlett(tmp.OLS$residuals)
    kpss.y <- KPSS(tmp.OLS$residuals, lr.var.y)

    DU1 <- c(rep(0, tb1), rep(1, n.obs - tb1))

    x <- cbind(
        rep(1, n.obs - 1),
        if (const) diff(DU1) else NULL,
        DU1[2:n.obs]
    )

    tmp.OLS <- OLS(d.y, x)

    lr.var.dy <- lr.var.bartlett(tmp.OLS$residuals)
    kpss.dy <- KPSS(tmp.OLS$residuals, lr.var.dy)

    lambda.kpss <- exp(-((500 * kpss.y * kpss.dy)^2))
    t.lambda.kpss <- lambda.kpss * t0 + m.ksi * (1 - lambda.kpss) * t1

    return(t.lambda.kpss)
}
