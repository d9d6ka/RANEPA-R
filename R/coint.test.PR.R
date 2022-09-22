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

    y_d <- GLS(y, z, opt.cbar)$residuals
    x_d <- GLS(x, z, opt.cbar)$residuals

    model <- OLS(y_d, x_d)
    ud.hat <- cbind(model$residuals)

    result <- resid.tests.PR(ud.hat, kmin, kmax, opt.cbar, det.comp)

    return(result)
}


resid.tests.PR <- function(ud.hat, kmin, kmax, opt.cbar, det.comp) {
    result <- list()

    n.obs <- nrow(ud.hat)

    if (ncol(ud.hat) > 1) {
        stop("ERROR:")
    }

    gls.tests <- matrix(0, 7, 1)
    min.bic <- Inf

    ud.hat.1 <- cbind(ud.hat[(1:(nrow(ud.hat) - 1)), 1])

    model.1 <- OLS(
        as.matrix(ud.hat[(2:nrow(ud.hat)), 1]),
        ud.hat.1
    )

    alpha.1 <- cbind(model.1$beta)
    ee.1 <- cbind(model.1$residuals)

    aa <- solve(t(ud.hat.1) %*% ud.hat.1)
    ss <- c(t(ee.1) %*% ee.1) / (nrow(ee.1) - ncol(ud.hat.1))
    t.alpha <- (alpha.1 - 1) / sqrt(ss * aa[1, 1])
    sum.ud.hat <- sum(ud.hat[(1:(n.obs - 1)), 1]^2) / (n.obs - 1)^2

    lag.bic <- kmin
    while (lag.bic <= kmax) {
        dud.hat <- diffn(ud.hat, na = 0)
        reg <- lagn(ud.hat, 1, na = 0)

        h <- 1
        while (h <= lag.bic) {
            reg <- cbind(reg, lagn(dud.hat, h, na = 0))
            h <- h + 1
        }

        dud.hat.1 <- as.matrix(dud.hat[(lag.bic + 2):nrow(dud.hat), ])
        reg.1 <- as.matrix(reg[(lag.bic + 2):nrow(reg), ])

        model.2 <- OLS(dud.hat.1, reg.1)

        ee.2 <- cbind(model.2$residuals)
        ss.2 <- c(t(ee.2) %*% ee.2) / (nrow(ee.2) - ncol(reg.1))
        xx.2 <- solve(t(reg.1) %*% reg.1)

        if (lag.bic == 0) {
            sumb <- 0
        } else {
            sumb <- sum(model.2$beta[2:(lag.bic + 1)])
        }

        ss.adj <- ss.2 / ((1 - sumb)^2)

        cur.bic <- log(c(t(ee.2) %*% ee.2) / (n.obs - kmax)) +
            (log(n.obs - kmax) * (lag.bic)) / (n.obs - kmax)
        if (cur.bic < min.bic) {
            min.bic <- cur.bic
            kbic <- lag.bic

            gls.tests[1, 1] <- ((ud.hat[n.obs, 1]^2 / n.obs) - ss.adj) *
                ((2 * sum.ud.hat))^(-1)
            gls.tests[2, 1] <- sqrt(sum.ud.hat / ss.adj)
            gls.tests[3, 1] <- gls.tests[1, 1] * gls.tests[2, 1]
            gls.tests[4, 1] <- model.2$coefficients[1] / sqrt(ss.2 * xx.2[1, 1])
            gls.tests[5, 1] <- (n.obs - 1) * (c(alpha.1) - 1) -
                0.5 * (ss.adj - ss) * sum.ud.hat^(-1)
            gls.tests[6, 1] <- sqrt(ss / ss.adj) * t.alpha -
                (ss.adj - ss) / sqrt(4 * ss.adj * sum.ud.hat)
        }

        lag.bic <- lag.bic + 1
    }

    if (det.comp == 1 || det.comp == 3) {
        gls.tests[7, 1] <- ((opt.cbar^2) * (sum.ud.hat) -
            (opt.cbar * (ud.hat[n.obs, 1]^2 / n.obs))) / s22
    } else if (det.comp == 2) {
        gls.tests[7, 1] <- ((opt.cbar^2) * (sum.ud.hat) +
            (1 - opt.cbar) * (ud.hat[n.obs, 1]^2 / n.obs)) / s22
    } else {
        stop("ERROR:")
    }

    result$gls.tests <- gls.tests
    result$kbic <- kbic

    return(result)
}
