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
    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
    }

    N <- nrow(y)

    if (model < 0 & model > 6) {
        stop("ERROR: Try to specify the deterministic component again")
    }

    if (weakly.exog) {
        if (model == 0) {
            xt <- x
        } else if (1 <= model & model <= 4) {
            deter <- determinants.KPSS.1.break(model, N, break.point)
            xt <- cbind(deter, x)
        } else if (model == 5) {
            deter <- determinants.KPSS.1.break(1, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        } else if (model == 6) {
            deter <- determinants.KPSS.1.break(4, N, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        c(beta, resid, ., t.beta) %<-% OLS(y, xt)
    } else {
        bic.min <- Inf
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

    test <- KPSS(resid, lr.var.AK(resid))

    return(
        list(
            beta = beta,
            test = test,
            resid = resid,
            t.beta = t.beta,
            break.point = break.point
        )
    )
}


#' @title
#' KPSS-test of cointegration
#'
#' @description
#' Procedure for testing the null of cointegration in the possible presence of
#' structural breaks.
#'
#' @details
#' Computes the cointegration test with one unknown structural break
#' where the break point is estimated either minimizing the value of
#' the statistic or the sum of the squared residuals.
#' The estimation of the cointegrating relationship bases on DOLS.
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y (Tx1)-vector of the dependent variable
#' @param x (Txk)-matrix of explanatory stochastic regressors
#' @param model Scalar, denotes the deterministic component:
#' \describe{
#' \item{1}{for model An.}
#' \item{2}{for model A.}
#' \item{3}{for model B.}
#' \item{4}{for model C.}
#' \item{5}{for model D.}
#' \item{6}{for model E.}
#' }
#' @param weakly.exog Exogeneity of the stochastic regressors
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous
#' (DOLS is used in this case).}
#' }
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return (2x2)-matrix, where the first rows gives the value of
#' the min(SC) test and the estimated break point;
#' the second row gives the value of the SC statistic,
#' where the break point is estimated as min(SSR).
#'
#' @importFrom zeallot %<-%
#' @export
KPSS.1.break.unknown <- function(y, x, model, weakly.exog, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    temp.result <- matrix(data = 0, nrow = N - 5, ncol = 2)

    for (i in 3:(N - 3)) {
        if (ll.init + 2 < i & i < N - 5 - ll.init) {
            c(beta, tests, resid, t_b, tb) %<-%
                KPSS.1.break(y, x, model, i, weakly.exog, ll.init)
            temp.result[i - 2, 1] <- tests
            temp.result[i - 2, 2] <- drop(t(resid) %*% resid)
        } else {
            temp.result[i - 2, 1] <- 2^20
            temp.result[i - 2, 2] <- 2^20
        }
    }

    minSC <- min(temp.result[, 1])
    tbe <- which.min(temp.result[, 1])
    result <- cbind(minSC, 2 + tbe)

    tbe <- which.min(temp.result[, 2])
    minSC <- temp.result[tbe, 1]
    result <- rbind(
        result,
        cbind(minSC, 2 + tbe)
    )

    colnames(result) <- c("stat", "tb")
    return(result)
}
