#' @title
#' Calculating long-run variance
#'
#' @param e (Tx1) vector or residuals.
#'
#' @return Long-run variance.
lr.var.bartlett <- function(e, l = NULL) {
    if (!is.matrix(e)) e <- as.matrix(e)

    N <- nrow(e)

    if (is.null(l))
        l <- trunc(4 * ((N / 100)^(1 / 4)))

    # lrv <- drop(t(e) %*% e) / N
    acf.e <- ACF(e)
    lrv <- acf.e[1]
    for (i in 1:l) {
        W <- (1 - i / (l + 1))
        lrv <- lrv + 2 * W * acf.e[1 + i]
    }
    return(lrv)
}


#' Estimating heteroscedasticity and autocorrelation consistent variance
#'
#' @references
#' Andrews, Donald W. K.
#' “Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation.”
#' Econometrica 59, no. 3 (1991): 817–58.
#' https://doi.org/10.2307/2938229.
#'
#' @importFrom zeallot %<-%
lr.var.quadratic <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    c(a, ., ., .) %<-% OLS(y[2:N, ], y[1:(N-1), ])
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
#' @param e (Tx1) vector or residuals.
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
#'
#' @return Long-run variance.
lr.var.bartlett.AK <- function(e) {
    if (!is.matrix(e)) e <- as.matrix(e)

    N <- nrow(e)
    k <- 0.8
    c(rho, ., ., .) %<-%
        OLS(
            e[2:N, , drop = FALSE],
            e[1:(N - 1), , drop = FALSE]
        )
    a <- drop(rho)
    l <- min(
        1.1447 * (4 * a^2 * N / ((1 + a)^2 * (1 - a)^2))^(1 / 3),
        1.1447 * (4 * k^2 * N / ((1 + k)^2 * (1 - k)^2))^(1 / 3)
    )
    l <- trunc(l)
    # lrv <- drop(t(e) %*% e) / N
    acf.e <- ACF(e)
    lrv <- acf.e[1]
    for (i in 1:l) {
        W <- (1 - i / (l + 1))
        lrv <- lrv + 2 * W * acf.e[1 + i]
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
#' @param e (Tx1) vector or residuals.
#' @param max.lag Maximum number of lags.
#' The exact number is selected by information criterions.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#' @param criterion The information crietreion: bic, aic or lwz.
#'
#' @return Long-run variance.
lr.var.SPC <- function(e,
                       max.lag = 0,
                       kernel = "bartlett",
                       criterion = "bic") {
    if (!is.matrix(e)) e <- as.matrix(e)
    if (max.lag < 0) max.lag <- 0
    if (!kernel %in% c("bartlett", "quadratic")) {
        warning("WARNING! Unknown kernel, Barlett is used")
        kernel <- "bartlett"
    }
    if (!criterion %in% c("bic", "aic", "lwz")) {
        warning("WARNING! Unknown criterion, BIC is used")
        criterion <- "bic"
    }

    N <- nrow(e)

    info.crit.min <- log(drop(t(e) %*% e) / (N - max.lag))

    c(rho, res, ., ., k) %<-%
        AR(e, NULL, max.lag, criterion)
    info.crit <- info.criterion(res, k)[[criterion]]

    if (info.crit.min < info.crit) {
        k <- 0
        rho <- 0
        res <- e
    }

    acf.res <- ACF(res)
    if (k == 0) {
        lrv <- acf.res[1]
        #lrv <- drop(t(res) %*% res) / N
    } else {
        temp <- cbind(res, lagn(res, 1))
        temp <- temp[2:nrow(res), , drop = FALSE]

        c(a, ., ., .) %<-%
            OLS(
                temp[, 1, drop = FALSE],
                temp[, 2, drop = FALSE]
            )
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
