#' @importFrom zeallot %<-%
#' @export
PY <- function(y,
               const = FALSE, trend = FALSE,
               criterion = "aic",
               trim = 0.15,
               max.lag) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (const && !trend) {
        VR <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.30, -4.39, -4.39, -4.34, -4.32,
            -4.45, -4.42, -4.33, -4.27, -4.27
        ))
        c.v <- matrix(c(
            1.60, 2.07, 3.33,
            1.52, 1.97, 3.24,
            1.41, 1.88, 3.05,
            1.26, 1.74, 3.12,
            0.91, 1.33, 2.83
        ), ncol = 3, byrow = TRUE)
    } else if (!const && trend) {
        VR <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.27, -4.41, -4.51, -4.55, -4.56,
            -4.57, -4.51, -4.38, -4.26, -4.26
        ))
        c.v <- matrix(c(
            1.52, 2.02, 3.37,
            1.40, 1.93, 3.27,
            1.28, 1.86, 3.20,
            1.13, 1.67, 3.06,
            0.74, 1.28, 2.61
        ), ncol = 3, byrow = TRUE)
    } else if (const && trend) {
        VR <- matrix(c(
            0, 1, 0, 0,
            0, 0, 0, 1
        ), nrow = 2, ncol = 4, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.38, -4.65, -4.78, -4.81, -4.90,
            -4.88, -4.75, -4.70, -4.41, -4.41
        ))
        c.v <- matrix(c(
            2.96, 3.55, 5.02,
            2.82, 3.36, 4.78,
            2.65, 3.16, 4.59,
            2.48, 3.12, 4.47,
            2.15, 2.79, 4.57
        ), ncol = 3, byrow = TRUE)
    } else {
        stop("ERROR! Unknown model")
    }

    N <- nrow(y)

    vect1 <- matrix(0, nrow = trunc((1 - 2 * trim) * N) + 2, ncol = 1)

    for (tb in (max(trunc(trim * N), max.lag + 2)):(trunc((1 - trim) * N))) {
        lambda <- tb / N

        DU <- c(rep(0, tb), rep(1, N - tb))
        DT <- DU * (1:N - tb)

        x <- cbind(
            rep(1, N),
            if (const) DU else NULL,
            1:N,
            if (trend) DT else NULL
        )

        k.hat <- max(1, lag.selection(y, x, criterion, max.lag))
        c(., resid, ., .) %<-% OLS(y, x)

        d.resid <- c(NA, diff(resid))

        x.resid <- lagn(resid, 1)
        for (l in 1:(f.hat - 1)) {
            x.resid <- cbind(x.resid, lagn(d.resid, l))
        }

        c(beta, e.hat, ., .) %<-% OLS(resid[k.hat:N, , drop = FALSE],
                                      x.resid[k.hat:N, , drop = FALSE])

        VCV <- eq.solve(t(x.resid) %*% x.resid) *
            drop(t(e.hat) %*% e.hat) / nrow(e.hat)

        a.hat <- beta[1]
        v.a.hat <- VCV[1, 1]
        tau <- (a.hat - 1) / sqrt(v.a.hat)

        #TODO: Доделать






    }


    return(
        list(
            wald = 0,
            critical.value = 0
        )
    )
}
