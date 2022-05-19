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
#' @param weakly.exog Exogeneity of the stochastic regressors
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous
#' (DOLS is used in this case).}
#' }
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinKPSS-test).}
#' \item{resid}{Residuals of the model.}
#' \item{t.beta}{Individual significance t-statistics.}
#' \item{break_point}{Break points.}
#' }
#'
#' @import MASS
#' @importFrom zeallot %<-%
#' @export
KPSS.1.break <- function(y, x,
                          model, break.point,
                          weakly.exog = TRUE, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x))
        if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    if (model < 0 & model > 6)
        stop("ERROR: Try to specify the deterministic component again")

    if (weakly.exog) {
        if (model == 0)
            xt <- x
        else if (1 <= model & model <= 4) {
            deter <- determinants.KPSS.1.break(model, N, break.point)
            xt <- cbind(deter, x)
        }
        else if (model == 5) {
            deter <- determinants.KPSS.1.break(1, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }
        else if (model == 6) {
            deter <- determinants.KPSS.1.break(4, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        c(beta, resid, ., t.beta) %<-% OLS(y, xt)
    }
    else {
        bic.min <- 100000000
        for (i in ll.init:1) {
            c(beta, resid, bic, t.beta) %<-%
                DOLS(y, x, model, break.point, i, i)
            if (bic < bic.min) {
                bic.min <- bic
                beta.min <- beta
                t.beta.min <- t.beta
                resid.min <- resid
            }
        }
        resid <- resid.min
        beta <- beta.min
        t.beta <- t.beta.min
    }

    S.t <- apply(resid, 2, cumsum)
    test <- N^(-2) * drop(t(S.t) %*% S.t) / alrvr(resid)

    return(
        list(
            beta   = beta,
            test   = test,
            resid  = resid,
            t.beta = t.beta,
            break.point = break.point
        )
    )
}
