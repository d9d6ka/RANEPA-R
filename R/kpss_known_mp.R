#' KPSS-test with 2 known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known
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
#' @param break_point Array of structural breaks.
#' @param trend Include trend if `TRUE`.
#' @param weakly.exog Boolean where we specify whether the stochastic regressors are exogenous or not
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous (DOLS is used in this case).}
#' }
#' @param kinit Scalar, defines the initial number of leads and lags for DOLS.
#' @param kmax scalar, with the maximum order of the parametric correction. The final order of the parametric correction is selected using the BIC information criterion.
#' @param kernel \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinkpss-test).}
#' \item{resid}{Residuals of the model.}
#' \item{break_point}{Break points.}
#' }
#'
#' @importFrom zeallot %<-%
#'
#' @export
kpss_known_mp <- function(y, x, model, break_point, const = FALSE, trend = FALSE, weakly.exog = TRUE, kinit, kmax, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    if (weakly.exog) {
        xt <- cbind(
            x,
            determi_kpss_mp(model, N, break_point, const, trend)
        )

        c(beta, resid, p, t_b) %<-% olsqr(y, xt)
    }
    else {
        bic_op <- 100000000
        for (i in kinit:1) {
            c(beta, resid, bic, t_b) %<-% dols_mp(y, x, model, break_point, const, trend, i, i)
            if (bic < bic_op) {
                bic_op <- bic
                beta_op <- beta
                t_b_op <- t_b
                resid_op <- resid
            }
        }
        resid <- resid_op
        beta <- beta_op
        t_b <- t_b_op
    }

    s_t <- apply(resid, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr_kernel(resid, kmax, kernel)
    else
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)

    return(
        list(
            beta = beta,
            test = test,
            resid = resid,
            break_point = break_point
        )
    )
}
