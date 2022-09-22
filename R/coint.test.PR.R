coint.test.PR <- function(y, x, det_comp, kmin = 0) {
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

    opt_cbar <- options[det_comp]

    kmax <- round(4 * (n.obs / 100)^0.25)

    z <- if (det_comp == 1 || det_comp == 3) {
        as.matrix(rep(1, n.obs))
    } else if (det_comp == 2) {
        cbind(rep(1, n.obs), (1:n.obs))
    }

    y_d <- GLS(y, z, opt_cbar)$residuals
    x_d <- GLS(x, z, opt_cbar)$residuals

    model <- lm.fit(x_d, y_d)
    ud_hat <- cbind(model$residuals)

    result <- resid.tests.PR(ud_hat, kmin, kmax, opt_cbar, det_comp)

    return(result)
}


resid.tests.PR <- function(ud_hat, kmin, kmax, opt_cbar, det_comp) {
    result <- list()

    n.obs <- nrow(ud_hat)

    if (ncol(ud_hat) > 1) {
        stop("ERROR:")
    }

    gls_tests <- matrix(0, 7, 1)
    minbicgls <- Inf

    lm_model <- lm.fit(
        as.matrix(ud_hat[(1:(nrow(ud_hat) - 1)), 1]),
        as.matrix(ud_hat[(2:nrow(ud_hat)), 1])
    )
    alpha_1 <- cbind(lm_model$coefficients)
    ee_1 <- cbind(lm_model$residuals)

    ud_hat_1 <- cbind(ud_hat[(1:(nrow(ud_hat) - 1)), 1])

    aa <- solve(t(ud_hat_1) %*% ud_hat_1)
    ss <- c(t(ee_1) %*% ee_1) / (nrow(ee_1) - ncol(ud_hat_1))
    talpha <- (alpha_1 - 1) / sqrt(ss * aa[1, 1])
    sum_ud_hat <- sum(ud_hat[(1:(n.obs - 1)), 1]^2) / (n.obs - 1)^2

    lag_bic <- kmin
    while (lag_bic <= kmax) {
        dud_hat <- diffn(ud_hat, na = 0)
        reg <- lagn(ud_hat, 1, na = 0)

        h <- 1
        while (h <= lag_bic) {
            reg <- cbind(reg, lagn(dud_hat, h, na = 0))
            h <- h + 1
        }

        dud_hat_1 <- as.matrix(dud_hat[(lag_bic + 2):nrow(dud_hat), ])
        reg_1 <- as.matrix(reg[(lag_bic + 2):nrow(reg), ])

        model_2 <- lm.fit(reg_1, dud_hat_1)

        e2 <- cbind(model_2$residuals)
        nef2 <- nrow(e2)
        s2e2 <- c(t(e2) %*% e2) / (nef2 - ncol(reg_1))
        xx2 <- solve(t(reg_1) %*% reg_1)

        if (lag_bic == 0) {
            sumb <- 0
        } else {
            sumb <- sum(model_2$coefficients[2:(lag_bic + 1)])
        }

        s22 <- s2e2 / ((1 - sumb)^2)

        bicgls <- log(c(t(e2) %*% e2) / (n.obs - kmax)) +
            (log(n.obs - kmax) * (lag_bic)) / (n.obs - kmax)
        if (bicgls < minbicgls) {
            minbicgls <- bicgls
            kbic <- lag_bic

            gls_tests[1, 1] <- ((ud_hat[n.obs, 1]^2 / n.obs) - s22) *
                ((2 * sum_ud_hat))^(-1)
            gls_tests[2, 1] <- sqrt(sum_ud_hat / s22)
            gls_tests[3, 1] <- gls_tests[1, 1] * gls_tests[2, 1]
            gls_tests[4, 1] <- model_2$coefficients[1] / sqrt(s2e2 * xx2[1, 1])
            gls_tests[5, 1] <- (n.obs - 1) * (c(alpha_1) - 1) -
                0.5 * (s22 - ss) * sum_ud_hat^(-1)
            gls_tests[6, 1] <- sqrt(ss / s22) * talpha -
                (s22 - ss) / sqrt(4 * s22 * sum_ud_hat)
        }

        lag_bic <- lag_bic + 1
    }

    if (det_comp == 1 || det_comp == 3) {
        gls_tests[7, 1] <- ((opt_cbar^2) * (sum_ud_hat) -
            (opt_cbar * (ud_hat[n.obs, 1]^2 / n.obs))) / s22
    } else if (det_comp == 2) {
        gls_tests[7, 1] <- ((opt_cbar^2) * (sum_ud_hat) +
            (1 - opt_cbar) * (ud_hat[n.obs, 1]^2 / n.obs)) / s22
    } else {
        stop("ERROR:")
    }

    result$gls_tests <- gls_tests
    result$kbic <- kbic

    return(result)
}
