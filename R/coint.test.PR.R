coint.test.PR <- function(y, x, det.comp, kmin = 0) {
    n.obs <- nrow(y)
    n.var <- ncol(x)

    if (n.var == 1) {
        options <- c(-13.75, -20.50, -13.50)
    } else if (n.var == 2) {
        options <- c(-18.25, -23.75, -18.00)
    } else if (n.var == 3) {
        options <- c(-22.25, -27.25, -23.00)
    } else if (n.var == 4) {
        options <- c(-26.25, -30.75, -26.00)
    } else {
        options <- c(-30.00, -33.75, -29.75)
    }

    opt.cbar <- options[det.comp]

    kmax <- round(4 * (n.obs / 100)^(1 / 4))

    z <- if (det.comp == 1 || det.comp == 3) {
        as.matrix(rep(1, n.obs))
    } else if (det.comp == 2) {
        cbind(rep(1, n.obs), (1:n.obs))
    }

    y.d <- GLS(y, z, opt.cbar)$residuals
    x.d <- GLS(x, z, opt.cbar)$residuals

    model <- OLS(y.d, x.d)
    ud.hat <- cbind(model$residuals)

    result <- resid.tests.PR(ud.hat, kmin, kmax, opt.cbar, det.comp)

    return(result)
}


resid.tests.PR <- function(ud, kmin, kmax, c.bar, det.comp) {
    n.obs <- nrow(ud)

    if (ncol(ud) > 1) {
        stop("ERROR:")
    }

    gls.tests <- matrix(0, 7, 1)

    lag.ud <- ud[1:(n.obs - 1), 1, drop = FALSE]
    d.ud <- diffn(ud, na = 0)
    sum.ud.sq <- drop(t(lag.ud) %*% lag.ud)

    model.1 <- OLS(
        ud[2:n.obs, 1, drop = FALSE],
        ud[1:(n.obs - 1), 1, drop = FALSE]
    )

    rho.hat <- as.matrix(model.1$beta)
    omega <- as.matrix(model.1$residuals)
    s2.ud <- c(t(omega) %*% omega) / (nrow(omega) - 1)
    t.rho <- (rho.hat - 1) / sqrt(s2.ud / sum.ud.sq)

    min.bic <- Inf
    min.lag <- kmin
    lag.bic <- kmin
    while (lag.bic <= kmax) {
        tmp.reg <- lagn(ud, 1, na = 0)

        h <- 1
        while (h <= lag.bic) {
            tmp.reg <- cbind(
                tmp.reg,
                lagn(d.ud, h, na = 0)
            )
            h <- h + 1
        }

        tmp.reg <-
            tmp.reg[(lag.bic + 2):nrow(tmp.reg), , drop = FALSE]

        model.2 <- OLS(
            d.ud[(lag.bic + 2):nrow(d.ud), , drop = FALSE],
            tmp.reg
        )

        eta <- cbind(model.2$residuals)
        s2.eta <- c(t(eta) %*% eta) / (nrow(eta) - ncol(tmp.reg))
        xtx.inv <- solve(t(tmp.reg) %*% tmp.reg)

        if (lag.bic == 0) {
            sumb <- 0
        } else {
            sumb <- sum(model.2$beta[2:(lag.bic + 1)])
        }

        ss.adj <- s2.eta / ((1 - sumb)^2)

        cur.bic <- log(c(t(eta) %*% eta) / (n.obs - kmax)) +
            log(n.obs - kmax) * lag.bic / (n.obs - kmax)
        if (cur.bic < min.bic) {
            min.bic <- cur.bic
            min.lag <- lag.bic

            gls.tests[1, 1] <- (ud[n.obs, 1]^2 / n.obs - ss.adj) /
                (2 * sum.ud.sq / n.obs^2)
            gls.tests[2, 1] <- sqrt(2 * sum.ud.sq / (n.obs^2 * ss.adj))
            gls.tests[3, 1] <- gls.tests[1, 1] * gls.tests[2, 1]
            gls.tests[4, 1] <- model.2$beta[1] / sqrt(s2.eta * xtx.inv[1, 1])
            gls.tests[5, 1] <- (n.obs - 1) * (rho.hat - 1) -
                (ss.adj - s2.ud) / (2 * sum.ud.sq / n.obs^2)
            gls.tests[6, 1] <- sqrt(s2.ud / ss.adj) * t.rho -
                (ss.adj - s2.ud) / sqrt(4 * ss.adj * sum.ud.sq / n.obs^2)
        }

        lag.bic <- lag.bic + 1
    }

    if (det.comp == 1 || det.comp == 3) {
        gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / n.obs^2 -
            c.bar * ud[n.obs, 1]^2 / n.obs) / ss.adj
    } else if (det.comp == 2) {
        gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / n.obs^2 +
            (1 - c.bar) * ud[n.obs, 1]^2 / n.obs) / ss.adj
    } else {
        stop("ERROR:")
    }

    rownames(gls.tests) <- c(
        "MZ(rho)",
        "MSB",
        "MZ(t.rho)",
        "ADF",
        "Z(rho)",
        "Z(t.rho)",
        if (det.comp == 1 || det.comp == 3) {
            "MP(T, demeaned)"
        } else if (det.comp == 2) {
            "MP(T, detrended)"
        }
    )

    result <- list()
    result$gls.tests <- gls.tests
    result$kbic <- min.lag

    return(result)
}
