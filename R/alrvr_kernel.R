## Procedure ALRVR to estimate the long-run variance as
## in Sul, Phillips and Choi (2003).
## Quadratic Spectral and Bartlett kernels.
#' @import MASS
alrvr_kernel <- function(e, kmax = 0, kernel = "bartlett", criterion = "bic") {
    if (!is.matrix(e)) e <- as.matrix(e)
    if (kmax < 0) kmax <- 0
    if (! kernel %in% c("bartlett", "quadratic")) {
        warning("WARNING! Unknown kernel, Barlett is used")
        kernel <- "bartlett"
    }
    if (! criterion %in% c("bic", "aic", "lwz")) {
        warning("WARNING! Unknown criterion, BIC is used")
        criterion <- "bic"
    }

    t <- nrow(e)

    min_bic <- log(as.numeric(t(e) %*% e) / (t - kmax))
    k <- 0
    rho <- 0
    res <- e


    for (i in 1:kmax) {
        if (kmax == 0) break
        
        temp <- e
        for (j in 1:i) temp <- cbind(temp, lagn(e, j))
        temp <- temp[(1 + i):nrow(temp), , drop = FALSE]
        
        x_tmp <- temp[, 2:ncol(temp), drop = FALSE]
        rho_temp <- qr.solve(t(x_tmp) %*% x_tmp) %*%
            t(x_tmp) %*% temp[, 1, drop = FALSE]
        res_temp <- temp[, 1, drop = FALSE] - x_tmp %*% rho_temp
        
        if (criterion == "bic") {
            bic <- log(as.numeric(t(res_temp) %*% res_temp) / (t - kmax)) +
                (i * log(t - kmax) / (t - kmax))
        } else if (criterion == "aic") {
            bic <- log(as.numeric(t(res_temp) %*% res_temp) / (t - kmax)) +
                2 * i / (t - kmax)
        } else if (criterion == "lwz") {
            bic <- log(as.numeric(t(res_temp) %*% res_temp) / (t - kmax)) +
                0.299 * i * (log(t - kmax))^2.1
        }

        if (bic < min_bic) {
            min_bic <- bic
            k <- i
            rho <- rho_temp
            res <- res_temp
        }
    }


    if (k == 0)
        lrv <- as.numeric(t(res) %*% res) / t
    else {
        temp <- cbind(res, lagn(res, 1))
        temp <- temp[2:nrow(res), , drop = FALSE]
        x_tmp <- temp[, 2, drop = FALSE]
        a <- as.numeric(
            qr.solve(t(x_tmp) %*% x_tmp) %*% t(x_tmp) %*% temp[, 1, drop = FALSE]
        )

        if (kernel == "bartlett")
            l <- 1.1447 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 3) # nolint
        else if (kernel == "quadratic")
            l <- 1.3221 * (4 * a^2 * t / ((1 + a)^2 * (1 - a)^2))^(1 / 5) # nolint
        l <- trunc(l)

        lrv <- as.numeric(t(res) %*% res) / t
        for (i in 1:l) {
            if (l == 0) break
            
            if (kernel == "bartlett")
                ## Bartlett kernel
                w <- (1 - i / (l + 1))
            else if (kernel == "quadratic")
                ## Quadratic spectral kernel
                w <- 25 / (12 * pi^2 * (i / l)^2) *
                    (sin(6 * pi * i / (l * 5)) / (6 * pi * i / (l * 5)) -
                        cos(6 * pi * i / (l * 5)))
            lrv <- lrv + 2 * as.numeric(t(res[1:(nrow(res) - i), , drop = FALSE]) %*%
                                        res[(1 + i):nrow(res), , drop = FALSE]) * w / t
        }
    }

    lrv_recolored <- lrv / (1 - sum(rho))^2

    lrv <- min(lrv_recolored, t * 0.15 * lrv)
    
    return(lrv)
}
