#' LM-test of cointegration
#'
#' @description
#' Procedure for testing the null of cointegration in the possible presence of structural breaks
#'
#' @details
#' Computes the cointegration test with one unknown structural break where the break point is estimated either minimizing the value of the statistic or the sum of the squared residuals. The estimation of the cointegrating relationship bases on DOLS.
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y (Tx1)-vector of the dependent variable
#' @param x (Txk)-matrix of explanatory stochastic regressors
#' @param model Scalar, denotes the deterministic component:
#' \describe{
#' \item{1}{for model An}
#' \item{2}{for model A}
#' \item{3}{for model B}
#' \item{4}{for model C}
#' \item{5}{for model D}
#' \item{6}{for model E}
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
#' @return vec_out (2x2)-matrix, where the first rows gives the value of the min(SC) test and the estimated break point; the second row gives the value of the SC statistic, where the break point is estimated as min(SSR).
#'
#' @importFrom zeallot %<-%
#' @export
coi_end <- function(y, x, model, k2, cri) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    m_SC <- matrix(data = 0, nrow = t - 5, ncol = 2)

    for (i in 3:(t-3)) {
        if (cri[2] + 2 < i & i < t - 5 - cri[2]) {
            c(beta, tests, u, t_b) %<-% coi_kpss(y, x, model, i, k2, cri)
            m_SC[i - 2, 1] <- tests
            m_SC[i - 2, 2] <- as.numeric(t(u) %*% u)
        }
        else {
            m_SC[i - 2, 1] <- 2 ^ 20
            m_SC[i - 2, 2] <- 2 ^ 20
        }
    }

    minSC <- min(m_SC[, 1])
    tbe <- which.min(m_SC[, 1])
    vec_out <- cbind(minSC, 2 + tbe)

    tbe <- which.min(m_SC[, 2])
    minSC <- m_SC[tbe, 1]
    vec_out <- rbind(
        vec_out,
        cbind(minSC, 2 + tbe)
    )

    colnames(vec_out) <- c("stat", "tb")
    return(vec_out)
}
