#' KPSS-test with 2 known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known.
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
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
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinKPSS-test).}
#' \item{resid}{Residuals of the model.}
#' \item{t.beta}{t-statistics for `beta`.}
#' \item{DOLS.lags}{The estimated number of lags and leads in DOLS.}
#' \item{break_point}{Break points.}
#' }
#'
#' @importFrom zeallot %<-%
#'
#' @export
KPSS.N.breaks <- function(y, x,
                          model, break.point,
                          const = FALSE, trend = FALSE,
                          weakly.exog = TRUE,
                          ll.init, corr.max, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    if (weakly.exog) {
        xt <- cbind(
            x,
            determinants.KPSS.N.breaks(model, N, break.point, const, trend)
        )

        c(beta, resid, ., t.beta) %<-% OLS(y, xt)
        DOLS.lags <- 0
    }
    else {
        bic.min <- 100000000
        for (i in ll.init:1) {
            c(beta, resid, bic, t.beta) %<-%
                DOLS.N.breaks(y, x, model, break.point, const, trend, i, i)
            if (bic < bic.min) {
                bic.min <- bic
                beta.min <- beta
                t.beta.min <- t.beta
                resid.min <- resid
                DOLS.lags <- i
            }
        }
        resid <- resid.min
        beta <- beta.min
        t.beta <- t.beta.min
    }

    S.t <- apply(resid, 2, cumsum)
    if (!is.null(kernel))
        test <- N^(-2) * drop(t(S.t) %*% S.t) / alrvr.kernel(resid, corr.max, kernel)
    else
        test <- N^(-2) * drop(t(S.t) %*% S.t) / alrvr(resid)

    return(
        list(
            beta   = beta,
            test   = test,
            resid  = resid,
            t.beta = t.beta,
            DOLS.lags = DOLS.lags,
            break.point = break.point
        )
    )
}
