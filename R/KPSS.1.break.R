#' @title
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
#' @param y An input (LHS) time series of interest.
#' @param x A matrix of (RHS) explanatory stochastic regressors.
#' @param model A scalar equal to
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param break.point Position of the break point.
#' @param weakly.exog Exogeneity of the stochastic regressors
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return A list of:
#' * `beta`: DOLS estimates of the coefficients regressors,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: Individual significance t-statistics,
#' * `break_point`: Break points.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
KPSS.1.break <- function(y, x,
                         model, break.point,
                         weakly.exog = TRUE, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
    }

    n.obs <- nrow(y)

    if (model < 0 && model > 6) {
        stop("ERROR: Try to specify the deterministic component again")
    }

    if (weakly.exog) {
        if (model == 0) {
            xt <- x
        } else if (1 <= model && model <= 4) {
            deter <- determinants.KPSS.1.break(model, n.obs, break.point)
            xt <- cbind(deter, x)
        } else if (model == 5) {
            deter <- determinants.KPSS.1.break(1, n.obs, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        } else if (model == 6) {
            deter <- determinants.KPSS.1.break(4, n.obs, break.point)
            xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
            xt <- cbind(deter, x, xdu)
        }

        tmp.OLS <- OLS(y, xt)
        beta <- tmp.OLS$beta
        resids <- tmp.OLS$residuals
        t.beta <- tmp.OLS$t.beta
        rm(tmp.OLS)
    } else {
        bic.min <- Inf
        for (i in ll.init:1) {
            tmp.DOLS <- DOLS(y, x, model, break.point, i, i)
            beta <- tmp.DOLS$beta
            resids <- tmp.DOLS$residuals
            bic <- tmp.DOLS$bic
            t.beta <- tmp.DOLS$t.beta
            rm(tmp.DOLS)

            if (bic < bic.min) {
                bic.min <- bic
                beta.min <- beta
                t.beta.min <- t.beta
                resid.min <- resids
            }
        }
        resids <- resid.min
        beta <- beta.min
        t.beta <- t.beta.min
    }

    test <- KPSS(resids, lr.var.bartlett.AK(resids))

    return(
        list(
            beta = beta,
            test = test,
            residuals = resids,
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
#' @param model A scalar equal to
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param weakly.exog Exogeneity of the stochastic regressors
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return (2x2)-matrix, where the first rows gives the value of
#' the min(SC) test and the estimated break point;
#' the second row gives the value of the SC statistic,
#' where the break point is estimated as min(SSR).
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
KPSS.1.break.unknown <- function(y, x, model, weakly.exog, ll.init) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    n.obs <- nrow(y)

    temp.result <- matrix(data = 0, nrow = n.obs - 5, ncol = 2)

    for (i in 3:(n.obs - 3)) {
        if (ll.init + 2 < i && i < n.obs - 5 - ll.init) {
            tmp.kpss <- KPSS.1.break(y, x, model, i, weakly.exog, ll.init)
            temp.result[i - 2, 1] <- tmp.kpss$test
            temp.result[i - 2, 2] <-
                drop(t(tmp.kpss$residuals) %*% tmp.kpss$residuals)
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
