#' @title
#' Calculating long-run variance with Bartlett kernel
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A series of interest.
#' @param l Number of lags/leads. If `NULL` then it's estimated.
#'
#' @return Long-run variance.
lr.var.bartlett <- function(y, l = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    if (is.null(l)) {
        l <- trunc(4 * ((N / 100)^(1 / 4)))
    }

    # lrv <- drop(t(e) %*% e) / N
    acf.y <- ACF(y)
    lrv <- acf.y[1]
    for (i in 1:l) {
        W <- (1 - i / (l + 1))
        lrv <- lrv + 2 * W * acf.y[1 + i]
    }
    return(lrv)
}


#' @title
#' Calculating long-run variance with quadratic kernel
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A series of interest.
#'
#' @references
#' Andrews, Donald W. K.
#' “Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation.”
#' Econometrica 59, no. 3 (1991): 817–58.
#' https://doi.org/10.2307/2938229.
lr.var.quadratic <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    a <- OLS(y[2:N, ], y[1:(N-1), ])$beta
    a <- drop(a)

    acf.y <- ACF(y)

    lambda <- matrix(0, nrow = N - 1, ncol = 1)
    s <- as.matrix(1:(N - 1))

    m <- 1.3221 * (4 * a^2 * N / ((1 - a)^4))^(1 / 5)

    delta <- (6 * pi * s) / (5 * m)

    for (i in 1:(N - 1)) {
        lambda[i] <- (3 / delta[i]^2) *
            (sin(delta[i]) / delta[i] - cos(delta[i]))
    }

    return(
        list(
            lrv = drop(acf.y[1] + 2 * t(lambda) %*% acf.y[2:N]),
            m = m
        )
    )
}


#' @title
#' Calculating long-run variance
#'
#' @description
#' Procedure ALRVR to estimate the long-run variance
#' as in Andrews (1991) and Kurozumi (2002).
#'
#' @param y (Tx1) vector or residuals.
#'
#' @return Long-run variance.
#'
#' @references
#' Andrews, Donald W. K.
#' “Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation.”
#' Econometrica 59, no. 3 (1991): 817–58.
#' https://doi.org/10.2307/2938229.
#'
#' Kurozumi, Eiji.
#' “Testing for Stationarity with a Break.”
#' Journal of Econometrics 108, no. 1 (May 1, 2002): 63–99.
#' https://doi.org/10.1016/S0304-4076(01)00106-3.
lr.var.bartlett.AK <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)
    k <- 0.8
    rho <- OLS(
        y[2:N, , drop = FALSE],
        y[1:(N - 1), , drop = FALSE]
    )$beta
    a <- drop(rho)
    l <- min(
        1.1447 * (4 * a^2 * N / ((1 + a)^2 * (1 - a)^2))^(1 / 3),
        1.1447 * (4 * k^2 * N / ((1 + k)^2 * (1 - k)^2))^(1 / 3)
    )
    l <- trunc(l)
    # lrv <- drop(t(e) %*% e) / N
    acf.y <- ACF(y)
    lrv <- acf.y[1]
    for (i in 1:l) {
        W <- (1 - i / (l + 1))
        lrv <- lrv + 2 * W * acf.y[1 + i]
    }
    return(lrv)
}


#' @title
#' Calculating long-run variance
#'
#' @description
#' Procedure ALRVR to estimate the long-run variance as
#' in Sul, Phillips and Choi (2003).
#'
#' @details Used are Quadratic Spectral and Bartlett kernels.
#'
#' @param y (Tx1) vector or residuals.
#' @param max.lag Maximum number of lags.
#' The exact number is selected by information criterions.
#' @param kernel Kernel for calculating long-run variance
#' * `"bartlett""`: for Bartlett kernel,
#' * `"quadratic"`: for Quadratic Spectral kernel,
#' * `"NULL"`: for the Kurozumi's proposal, using Bartlett kernel.
#' @param criterion The information crietreion: bic, aic or lwz.
#'
#' @return Long-run variance.
#'
#' @references
#' Sul, Donggyu, Peter C. B. Phillips, and Chi-Young Choi.
#' “Prewhitening Bias in HAC Estimation.”
#' Oxford Bulletin of Economics and Statistics 67, no. 4 (August 2005): 517–46.
#' https://doi.org/10.1111/j.1468-0084.2005.00130.x.
lr.var.SPC <- function(y,
                       max.lag = 0,
                       kernel = "bartlett",
                       criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (max.lag < 0) max.lag <- 0
    if (!kernel %in% c("bartlett", "quadratic")) {
        warning("WARNING! Unknown kernel, Barlett is used")
        kernel <- "bartlett"
    }
    if (!criterion %in% c("bic", "aic", "lwz")) {
        warning("WARNING! Unknown criterion, BIC is used")
        criterion <- "bic"
    }

    N <- nrow(y)

    info.crit.min <- log(drop(t(y) %*% y) / (N - max.lag))

    tmp.AR <- AR(y, NULL, max.lag, criterion)
    rho <- tmp.AR$beta
    res <- tmp.AR$residuals
    k <- tmp.AR$lag
    rm(tmp.AR)
    info.crit <- info.criterion(res, k)[[criterion]]

    if (info.crit.min < info.crit) {
        k <- 0
        rho <- 0
        res <- y
    }

    acf.res <- ACF(res)
    if (k == 0) {
        lrv <- acf.res[1]
        #lrv <- drop(t(res) %*% res) / N
    } else {
        temp <- cbind(res, lagn(res, 1))
        temp <- temp[2:nrow(res), , drop = FALSE]

        a <- OLS(
            temp[, 1, drop = FALSE],
            temp[, 2, drop = FALSE]
        )$beta
        a <- drop(a)

        if (kernel == "bartlett") {
            l <- 1.1447 * (4 * a^2 * N / ((1 + a)^2 * (1 - a)^2))^(1 / 3)
        } else if (kernel == "quadratic") {
            l <- 1.3221 * (4 * a^2 * N / ((1 - a)^4))^(1 / 5)
        }
        l <- trunc(l)

        lrv <- acf.res[1]
        for (i in 1:l) {
            if (l == 0) break

            if (kernel == "bartlett") {
                ## Bartlett kernel
                W <- (1 - i / (l + 1))
            } else if (kernel == "quadratic") {
                ## Quadratic spectral kernel
                delta <- 6 * pi * (i / l) / 5
                W <- (3 / delta^2) * (sin(delta) / delta - cos(delta))
            }
            lrv <- lrv + 2 * W * acf.res[1 + i]
        }
    }

    lrv.recolored <- lrv / (1 - sum(rho))^2

    lrv <- min(lrv.recolored, N * 0.15 * lrv)

    return(lrv)
}
