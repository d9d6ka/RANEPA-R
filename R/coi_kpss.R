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
#' @param y (Tx1)-vector of the dependent variable
#' @param x (Txk)-matrix of explanatory stochastic regressors
#' @param model (2x1)-vector, denotes the deterministic component:
#' \describe{
#' \item{model[1]}{\describe{
#' \item{1}{for model An}
#' \item{2}{for model A}
#' \item{3}{for model B}
#' \item{4}{for model C}
#' \item{5}{for model D}
#' \item{6}{for model E}}
#' }
#' \item{model[2]}{Position of the break point.}
#' }
#' @param k2 Scalar. Indicates the deterministic component of the stochastic regressors:
#' \describe{
#' \item{1}{if x is a matrix of I(1) variables with drift.}
#' \item{2}{if x is a matrix of I(1) variables with a quadratic trend.}
#' }
#' @param cri (2x1)-vector where we specify whether the stochastic regressors are exogenous or not, and the initial number of leads and lags for the DOLS estimation.
#' \describe{
#' \item{cri[1]}{\describe{
#' \item{0}{if the regressors are weakly exogenous,}
#' \item{1}{if the regressors are not weakly exogenous (DOLS is used in this case).}
#' }
#' }
#' \item{cri[2]}{Scalar, defines the initial number of leads and lags for DOLS}
#' }
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinkpss-test).}
#' \item{u}{Residuals of the model}
#' \item{t_b}{Individual significance t-statistics}
#' }
#'
#' @importFrom zeallot %<-%
#' @export
coi_kpss <- function(y, x, model, tb, k2, cri) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    if (model < 0 & model > 6)
        stop("ERROR: Try to specify the deterministic component again")

    if (cri[1] == 0) {
        if (model == 0)
            xt <- x
        else if (1 <= model & model <= 4) {
            deter <- determi(model, t, tb)
            xt <- cbind(deter, x)
        }
        else if (model == 5) {
            deter <- determi(1, t, tb)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }
        else if (model == 6) {
            deter <- determi(4, t, tb)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        beta <- solve(t(xt) %*% xt) %*% t(xt) %*% y
        u <- y - xt %*% beta
        t_b <- tb
    }
    else if (cri[1] == 1) {
        k <- cri[2]
        bic_op <- 100000000
        for (i in k:1) {
            c(beta, u, bic, t_b) %<-% dols(y, x, model, k, k, tb)
            if (bic < bic_op) {
                bic_op <- bic
                beta_op <- beta
                t_b_op <- t_b
                u_op <- u
            }
        }
        u <- u_op
        beta <- beta_op
        t_b <- t_b_op
    }

    sg <- alrvr(u)
    tests <- solve((t^2) * sg) * apply(apply(u, 2, cumsum)^2, 2, sum)

    return(
        list(
            beta = beta,
            tests = tests,
            u = u,
            t_b = t_b
        )
    )
}
