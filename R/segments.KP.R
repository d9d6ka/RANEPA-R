segments.KP <- function(y,
                        const = FALSE, trend = FALSE,
                        breaks = 1,
                        criterion = "aic",
                        trim = 0.15,
                        max.lag = 1) {
    if (const && !trend) {
        R <- matrix(c(0, 1, 0), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.30, -4.39, -4.39, -4.34, -4.32,
            -4.45, -4.42, -4.33, -4.27, -4.27
        ))
    } else if (!const && trend) {
        R <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.27, -4.41, -4.51, -4.55, -4.56,
            -4.57, -4.51, -4.38, -4.26, -4.26
        ))
    } else if (const && trend) {
        R <- matrix(c(
            0, 1, 0, 0,
            0, 0, 0, 1
        ), nrow = 2, ncol = 4, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.38, -4.65, -4.78, -4.81, -4.90,
            -4.88, -4.75, -4.70, -4.41, -4.41
        ))
    } else {
        stop("ERROR! Unknown model")
    }

    N <- nrow(y)
    x.const <- rep(1, N)
    x.trend <- 1:N

    h <- trunc(trim * N)

    if (l == 0) {
        date.vec <- c(1, N + 1)
    } else {
        SSR.data <- SSR.matrix(y, cbind(x.const, x.trend), h)
        dates <- segments.OLS(y, cbind(x.const, x.trend), l, h, SSR.data)
        date.vec <- c(1, drop(dates$break.point), N + 1)
    }

    res.wald <- -Inf

    for (i in 1:(l + 1)) {
        if vect1 <- rep(0, l + 1)
        t.low <- trunc(
            max(
                date.vec[i] + (date.vec[i + 1] - date.vec[i]) * trim,
                max.lag + 2
            )
        )
        t.high <-
            trunc(date.vec[i + 1] + (date.vec[i + 1] - date.vec[i]) * trim)
        if (t.low < t.high - 1) {
            for (tb in t.low:t.high) {
                DU <- c(rep(0, tb), rep(1, N - tb))
                DT <- DU * (1:N - tb)

                x <- cbind(
                    x.const,
                    x.trend - date.vec[i],
                    if (const) DU else NULL,
                    if (trend) DT else NULL
                )

                lambda <- tb / date.vec[i + 1]

                y.i <- y[(date.vec[i] + 1):date.vec[i + 1], , drop = FALSE]
                x.i <- x[(date.vec[i] + 1):date.vec[i + 1], , drop = FALSE]

                k.hat <- lag.selection(y.i, x.i, criterion, max.lag)

                c(., u, ., .) %<-% OLS(y.i, x.i)

                d.u <- as.matrix(c(0, diff(u)))
                y.u <- u[k.hat:N, , drop = FALSE]
                x.u <- lagn(u, 1, na = 0)
                for (j in 1:(k.hat - 1))
                    x.u <- cbind(x.u, lagn(d.u, j, na = 0))
                x.u <- x.u[k.hat:N, , drop = FALSE]

                c(., e.hat, ., .) %<-% OLS(y.u, x.u)

                VCV <- drop(t(e.hat) %*% e.hat) / nrow(e.hat) *
                    qr.solve(t(x.u) %*% x.u)

                a.hat <- beta[1]
                var.a.hat <- VCV[1, 1]
                tau <- (a.hat - 1) / sqrt(var.a.hat)

                tau05 <- v.t[ceiling(lambda * 10)]

                IP <- trunc((k.hat + 1) / 2)
                r <- ncol(x)
                k <- 10

                T.i <- date.vec[i + 1] - date.vec[i]
                c1 <- sqrt((1 + r) * T.i)
                c2 <- ((r + 1) * T.i - tau05^2 * (IP + T.i)) /
                    (tau05 * (tau05 + k) * (IP + T.i))

                if (tau > tau05)
                    c.tau <- -tau
                if (tau <= tau05 && tau > -k)
                    c.tau <- IP * tau / N - (r + 1) / (tau + c2 * (tau + k))
                if (tau <= -k && tau > -c1)
                    c.tau <- IP * tau / N - (r + 1) / tau
                if (tau <= -c1)
                    c.tau <- 0

                a.hat.M <- a.hat + c.tau * sqrt(var.a.hat)
                if (a.hat.M >= 1)
                    a.hat.M <- 1
                else if (abs(a.hat.M) < 1)
                    a.hat.M <- a.hat.M
                else
                    a.hat.M <- -0.99

                CR <- sqrt(N) * abs(a.hat.M - 1)
                if (CR <= 1) a.hat.M <- 1

            }
        }
    }

    return(0)
}
