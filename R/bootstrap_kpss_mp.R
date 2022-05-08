#' Calculating p-values using bootstrap
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known.
#'
#' See Cavaliere and Taylor (2006) for further details.
#'
#' @param y (Tx1)-vector of time series.
#' @param x (Txk)-matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param break.point Array of structural breaks.
#' @param trend Include trend if `TRUE`.
#' @param weakly.exog Boolean where we specify whether the stochastic regressors are exogenous or not
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous (DOLS is used in this case).}
#' }
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#' @param corr.max scalar, with the maximum order of the parametric correction. The final order of the parametric correction is selected using the BIC information criterion.
#' @param kernel \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#' @param iter Number of bootstrap iterations.
#' @param bootstrap Type of bootstrapping:
#' \describe{
#' \item{sample}{sampling from residuals with replacement.}
#' \item{Cavaliere-Taylor}{multiplying residuals by N(0, 1)-distributed variable.}
#' \item{Rademacher}{multiplying residuals by Rademacher-distributed variable.}
#' }
#'
#' @return List of 3 elements:
#' \describe{
#' \item{test}{The value of KPSS test statistic.}
#' \item{p_value}{The estimates p-value.}
#' \item{bootstrapped}{Bootstrapped auxiliary statistics.}}
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
bootstrap_kpss_mp <- function(y, x,
                              model, break.point,
                              const = FALSE, trend = FALSE,
                              weakly.exog = TRUE,
                              ll.init, corr.max, kernel,
                              iter = 9999,
                              bootstrap = "sample") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    c(., test, u, ., dols_lags, .) %<-%
        kpss_known_mp(
            y, x,
            model, break.point,
            const, trend,
            weakly.exog,
            ll.init, corr.max, kernel
        )

    if (weakly.exog) {
        xreg <- cbind(
            x,
            determi_kpss_mp(model, N, break.point, const, trend)
        )
    }
    else {
        c(., xreg) %<-%
            dols_prepare_mp(
                y, x,
                model, break.point,
                const, trend,
                dols_lags, dols_lags
            )
    }

    cores <- detectCores() # nolint

    progress_bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress_bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    registerDoSNOW(cluster)

    result <- foreach(
        i = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        if (bootstrap == "sample") {
            tmp_y <- sample(u, length(u), replace = TRUE)
        }
        else if (bootstrap == "Cavaliere-Taylor") {
            z <- rnorm(length(u))
            tmp_y <- z * u
        }
        else if (bootstrap == "Rademacher") {
            z <- sample(c(-1, 1), length(u), replace = TRUE)
            tmp_y <- z * u
        }

        c(beta, resid, ., t_beta) %<-% olsqr(tmp_y, xreg)

        s_t <- apply(resid, 2, cumsum)
        if (!is.null(kernel))
            tmp_test <- N^(-2) * drop(t(s_t) %*% s_t) /
                alrvr_kernel(resid, corr.max, kernel)
        else
            tmp_test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)
        tmp_test
    }

    stopCluster(cluster)

    p_value <- (1 / iter) * sum(I(test <= result))

    return(
        list(
            test = test,
            p_value = p_value,
            bootstrapped = result
        )
    )
}
