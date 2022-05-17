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
kpss.1.break.unknown <- function(y, x, model, weakly.exog, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    temp.result <- matrix(data = 0, nrow = N - 5, ncol = 2)

    for (i in 3:(N-3)) {
        if (ll.init + 2 < i & i < N - 5 - ll.init) {
            c(beta, tests, resid, t_b, tb) %<-%
                kpss.1.break(y, x, model, i, weakly.exog, ll.init)
            temp.result[i - 2, 1] <- tests
            temp.result[i - 2, 2] <- drop(t(resid) %*% resid)
        }
        else {
            temp.result[i - 2, 1] <- 2 ^ 20
            temp.result[i - 2, 2] <- 2 ^ 20
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
