#' @title
#' Sequential Perron-Yabu (2009) statistic for breaks at unknown date.
#'
#' @param y The input series of interest.
#' @param const Allowing the break in constant.
#' @param breaks Number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#' @param max.lag The maximum possible lag in the model.
#'
#' @import MASS
#' @importFrom zeallot %<-%
#'
#' @export
PY.sequential <- function(y,
                          const = FALSE,
                          breaks = 1,
                          criterion = "aic",
                          trim = 0.15,
                          max.lag = 1) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (!const) {
        R <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.27, -4.41, -4.51, -4.55, -4.56,
            -4.57, -4.51, -4.38, -4.26, -4.26
        ))
    } else if (const) {
        R <- matrix(c(
            0, 1, 0, 0,
            0, 0, 0, 1
        ), nrow = 2, ncol = 4, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.38, -4.65, -4.78, -4.81, -4.90,
            -4.88, -4.75, -4.70, -4.41, -4.41
        ))
    }

    N <- nrow(y)
    x.const <- rep(1, N)
    x.trend <- 1:N

    h <- trunc(trim * N)

    if (breaks == 0) {
        date.vec <- c(1, N + 1)
    } else {
        SSR.data <- SSR.matrix(y, cbind(x.const, x.trend), h)
        dates <- segments.OLS(y, cbind(x.const, x.trend), breaks, h, SSR.data)
        date.vec <- c(1, drop(dates$break.point), N + 1)
    }

    res.wald <- -Inf

    for (i in 1:(breaks + 1)) {
        vect1 <- rep(0, N)
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
                lambda <- tb / date.vec[i + 1]

                DU <- as.numeric(x.trend > tb)
                DT <- DU * (x.trend - tb)

                x <- cbind(
                    x.const,
                    x.trend - date.vec[i],
                    if (const) DU else NULL,
                    DT
                )

                y.i <- y[(date.vec[i] + 1):date.vec[i + 1], , drop = FALSE]
                x.i <- x[(date.vec[i] + 1):date.vec[i + 1], , drop = FALSE]

                k.hat <- max(1, lag.selection(y.i, x.i, max.lag, criterion))

                c(., resid, ., .) %<-% OLS(y.i, x.i)

                d.resid <- as.matrix(c(0, diff(resid)))

                y.u <- resid[k.hat:nrow(resid), , drop = FALSE]
                x.u <- lagn(resid, 1, na = 0)
                if (k.hat > 1) {
                    for (l in 1:(k.hat - 1)) {
                        x.u <- cbind(x.u, lagn(d.resid, l, na = 0))
                    }
                }
                x.u <- x.u[k.hat:N, , drop = FALSE]

                c(beta.u, u.resid, ., .) %<-% OLS(y.u, x.u)

                VCV <- qr.solve(t(x.u) %*% x.u) *
                    drop(t(u.resid) %*% u.resid) / nrow(u.resid)

                a.hat <- beta.u[1]
                var.a.hat <- VCV[1, 1]
                tau <- (a.hat - 1) / sqrt(var.a.hat)

                tau05 <- v.t[ceiling(lambda * 10)]

                IP <- trunc((k.hat + 1) / 2)

                k <- 10
                k.x <- ncol(x)

                T.i <- date.vec[i + 1] - date.vec[i]
                c1 <- sqrt((1 + k.i) * T.i)
                c2 <- ((1 + k.x) * T.i - tau05^2 * (IP + T.i)) /
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

                y.g <- rbind(
                    y[(1 + date.vec[i]), , drop = FALSE],
                    y[(2 + date.vec[i]):(date.vec[i + 1]), , drop = FALSE] -
                    a.hat.M * y[(1 + date.vec[i]):(date.vec[i + 1] - 1), ,
                                drop = FALSE]
                )
                x.g <- rbind(
                    x[(1 + date.vec[i]), , drop = FALSE],
                    x[(2 + date.vec[i]):(date.vec[i + 1]), , drop = FALSE] -
                    a.hat.M * x[(1 + date.vec[i]):(date.vec[i + 1] - 1), ,
                                drop = FALSE]
                )

                c(beta.g, g.resid, ., .) %<-% OLS(y.g, x.g)

                if (k.hat == 1) {
                    h0 <- drop(t(g.resid) %*% g.resid) / nrow(g.resid)
                } else {
                    if (a.hat.M == 1) {
                        x.v <- NULL
                        for (k.i in 1:(k.hat - 1))
                            x.v <- cbind(x.v, lagn(g.resid, k.i, na = 0))

                        y.v <- g.resid[(k.hat - 1):nrow(g.resid)]
                        x.v <- x.v[(k.hat - 1):nrow(g.resid), , drop = FALSE]

                        c(beta.v, v.resid, ., .) %<-% OLS(y.v, x.v)

                        if (!const) {
                            h0 <- (drop(t(v.resid) %*% v.resid) / (T.i - k.hat)) /
                                ((1 - sum(beta.v))^2)
                        }
                        if (const) {
                            BETAS <- matrix(0, nrow = k.hat - 1, ncol = 4)
                            for (k.i in 1:(k.hat - 1)) {
                                DU.ki <- as.numeric(x.trend > tb - k.i)
                                DT.ki <- DU.ki * (x.trend - tb)
                                x.ki <- cbind(
                                    x.const,
                                    DU.ki,
                                    x.trend,
                                    DT.ki
                                )
                                x.g.ki <- rbind(
                                    x.ki[(1 + date.vec[i]), ],
                                    x.ki[(2 + date.vec[i]):(date.vec[i + 1]), ] -
                                    a.hat.M * x.ki[(1 + date.vec[i]):(date.vec[i + 1] - 1), ]
                                )
                                c(beta.ki, ., ., .) %<-% OLS(y.g, x.g.ki)
                                BETAS[k.i, ] <- drop(beta.ki)
                                sig.e <- drop(t(v.resid) %*% v.resid) / (T.i - k.hat)
                                h0 <- sig.e / ((1 - sum(beta.v))^2)
                                beta.g[2] <- (sqrt(h0) / sqrt(sig.e)) *
                                    (beta.g[2] - drop(t(BETAS[, 2]) %*% beta.v))
                            }
                        }
                    }

                    if (abs(a.hat.M) < 1)
                        c(h0, m) <- h0W(g.resid)
                }

                VCV <- h0 * qr.solve(t(x.g) %*% x.g)
                vect1[tb] <- t(R %*% beta.g) %*%
                    qr.solve(R %*% VCV %*% t(R)) %*% (R %*% beta.g)
            }

            vect1 <- vect1[t.low:t:max]
            wald[i] <- log(sum(exp(vect1 / 2)) / T.i)
        }
    }

    return(
        list(
            wald = max(wald)
        )
    )
}
