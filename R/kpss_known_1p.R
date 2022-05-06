#' KPSS-test with known structural break
#'
#' @description
#' Computes the cointegration test with one known structural break.
#'
#' @details
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-matrix of explanatory stochastic regressors.
#' @param model \describe{
#' \item{1}{for model An.}
#' \item{2}{for model A.}
#' \item{3}{for model B.}
#' \item{4}{for model C.}
#' \item{5}{for model D.}
#' \item{6}{for model E.}
#' }
#' @param break.point Position of the break point.
#' @param weakly.exog Boolean where we specify whether the stochastic regressors are exogenous or not
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous (DOLS is used in this case).}
#' }
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinkpss-test).}
#' \item{resid}{Residuals of the model.}
#' \item{t_beta}{Individual significance t-statistics.}
#' \item{break_point}{Break points.}
#' }
#'
#' @import MASS
#' @importFrom zeallot %<-%
#' @export
kpss_known_1p <- function(y, x, model, break.point, weakly.exog = TRUE, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    if (model < 0 & model > 6)
        stop("ERROR: Try to specify the deterministic component again")

    if (weakly.exog) {
        if (model == 0)
            xt <- x
        else if (1 <= model & model <= 4) {
            deter <- determi_kpss_1p(model, N, break.point)
            xt <- cbind(deter, x)
        }
        else if (model == 5) {
            deter <- determi_kpss_1p(1, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }
        else if (model == 6) {
            deter <- determi_kpss_1p(4, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        c(beta, resid, p, t_b) %<-% olsqr(y, xt)
    }
    else {
        bic_min <- 100000000
        for (i in ll.init:1) {
            c(beta, resid, bic, t_beta) %<-% dols(y, x, model, i, i, break.point)
            if (bic < bic_min) {
                bic_min <- bic
                beta_min <- beta
                t_beta_min <- t_beta
                resid_min <- resid
            }
        }
        resid <- resid_min
        beta <- beta_min
        t_beta <- t_beta_min
    }

    s_t <- apply(resid, 2, cumsum)
    test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)

    return(
        list(
            beta   = beta,
            test   = test,
            resid  = resid,
            t_beta = t_beta,
            break.point = break.point
        )
    )
}
