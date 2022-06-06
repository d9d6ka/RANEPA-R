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
        if (k.hat > 1) {
            for (l in 1:(k.hat - 1)) {
                x.resid <- cbind(x.resid, lagn(d.resid, l))
            }
        }

        c(beta, e.hat, ., .) %<-% OLS(resid[(1 + k.hat):N, , drop = FALSE],
                                      x.resid[(1 + k.hat):N, , drop = FALSE])

        VCV <- qr.solve(t(x.resid[(1 + k.hat):N, , drop = FALSE]) %*% x.resid[(1 + k.hat):N, , drop = FALSE]) *
            drop(t(e.hat) %*% e.hat) / nrow(e.hat)

        a.hat <- beta[1]
        v.a.hat <- VCV[1, 1]
        tau <- (a.hat - 1) / sqrt(v.a.hat)

        tau.05 <- v.t[ceiling(lambda * 10)]

        IP <- trunc((k.hat + 1) / 2)
        k <- 10
        k.x <- ncol(x)

        c1 <- sqrt((1 + k.x) * N)
        c2 <- ((1 + k.x) * N - tau.05^2 * (IP + N)) /
            (tau.05 * (tau.05 + k.x) * (IP + N))

        if (tau > tau.05)
            c.tau <- -tau
        else if (tau <= tau.05 && tau > -k)
            c.tau <- IP * tau / N - (k.x + 1) / (tau + c2 * (tau + k))
        else if (tau <= -k && tau > -c1)
            c.tau <- IP * tau / N - (k.x + 1) / tau
        else if (tau <= -c1)
            c.tau <- 0

        rhomd1 <- a.hat + c.tau * sqrt(v.a.hat)
        if (rhomd1 >= 1)
            amu <- 1
        else if (abs(rhomd1) < 1)
            amu <- rhomd1
        else
            amu <- -0.99

        CR <- sqrt(N) * abs(amu - 1)
        if (CR <= 0) amu <- 1

        g.y <- rbind(y[1, , drop = FALSE], y[2:N, , drop = FALSE] - amu * y[1:(N - 1), , drop = FALSE])
        g.x <- rbind(x[1, , drop = FALSE], x[2:N, , drop = FALSE] - amu * x[1:(N - 1), , drop = FALSE])

        c(g.beta, v.resid, ., .) %<-% OLS(g.y, g.x)

        if (k.hat == 1)
            h0 <- drop(t(v.resid) %*% v.resid) / nrow(v.resid)
        else {
            if (amu == 1) {
                v.x <- NULL
                for (k.i in 1:(k.hat - 1)) v.x <- cbind(v.x, lagn(v.resid, k.i))
                v.y <- v.resid[k.hat:(N - 1), , drop = FALSE]
                c(beta, e.resid, ., .) %<-%
                    OLS(v.y,
                        v.x[k.hat:(N - 1), , drop = FALSE])

                if (const && !trend) {
                    v.beta <- matrix(0, nrow = k.hat - 1, ncol = 3)
                    for (k.i in 1:(k.hat - 1)) {
                        DU.ki <- c(rep(0, tb - k.i), rep(1, N - (tb - k.i)))
                        x.ki <- cbind(
                            rep(1, N),
                            1:N,
                            DU.ki
                        )
                        g.x.ki <- rbind(x.ki[1, ],
                                        x.ki[2:N, ] - amu * x.ki[1:(N - 1), ])
                        c(beta.ki, ., ., .) %<-% OLS(g.y, g.x.ki)
                        v.beta[k.i, ] <- drop(beta.ki)
                    }
                    g.beta[2] <- g.beta[2] - drop(t(v.beta[, 2, drop = FALSE]) %*% beta)
                    h0 <- drop(t(e.resid) %*% e.resid) / (N - k.hat)
                } else if (!const && trend) {
                    h0 <- (drop(t(e.resid) %*% e.resid) / (N - k.hat)) /
                        ((1 - sum(beta))^2)
                } else {
                    v.beta <- matrix(0, nrow = k.hat - 1, ncol = 4)
                    for (k.i in 1:(k.hat - 1)) {
                        DU.ki <- c(rep(0, tb - k.i), rep(1, N - (tb - k.i)))
                        DT.ki <- DU.ki * (1:N - tb)
                        x.ki <- cbind(
                            rep(1, N),
                            DU.ki,
                            1:N,
                            DT.ki
                        )
                        g.x.ki <- rbind(x.ki[1, ],
                                        x.ki[2:N, ] - amu * x.ki[1:(N - 1), ])
                        c(beta.ki, ., ., .) %<-% OLS(g.y, g.x.ki)
                        v.beta[k.i, ] <- drop(beta.ki)
                        sig.e <- drop(t(e.resid) %*% e.resid) / (N - k.hat)
                        h0 <- sige / ((1 - sum(beta))^2)
                        g.beta[2] <- (g.beta[2] - drop(t(v.beta[, 2, drop = FALSE]) %*% beta)) *
                            sqrt(h0) / sqrt(sig.e)
                    }
                }
            }

            if (abs(amu) < 1)
                c(h0, m) <- h0W(v.resid)
        }

        VCV = h0 * qr.solve(t(g.x) %*% g.x)
        vect1[tb - trunc(trim * N)] = t(VR %*% g.beta) %*%
            qr.solve(VR %*% VCV %*% t(VR)) %*% (VR %*% g.beta)
    }

    wald <- log(sum(exp(vect1 / 2)) / N)
    if (trim == 0.01) cv <- c.v[1, ]
	if (trim == 0.05) cv <- c.v[2, ]
	if (trim == 0.10) cv <- c.v[3, ]
	if (trim == 0.15) cv <- c.v[4, ]
	if (trim == 0.25) cv <- c.v[5, ]

    return(
        list(
            wald = wald,
            critical.value = cv
        )
    )
}
