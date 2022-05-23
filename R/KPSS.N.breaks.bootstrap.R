#' @title
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
#' @param const Include constant if **TRUE**.
#' @param trend Include trend if **TRUE**.
#' @param weakly.exog Exogeneity of the stochastic regressors
#' \describe{
#' \item{TRUE}{if the regressors are weakly exogenous,}
#' \item{FALSE}{if the regressors are not weakly exogenous
#' (DOLS is used in this case).}
#' }
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#' @param corr.max scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the BIC
#' information criterion.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#' @param iter Number of bootstrap iterations.
#' @param bootstrap Type of bootstrapping:
#' \describe{
#' \item{sample}{sampling from residuals with replacement.}
#' \item{Cavaliere-Taylor}{multiplying residuals by \eqn{N(0, 1)}-distributed
#' variable.}
#' \item{Rademacher}{multiplying residuals by Rademacher-distributed variable.}
#' }
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic or lwz.
#'
#' @return List of 3 elements:
#' \describe{
#' \item{test}{The value of KPSS test statistic.}
#' \item{p.value}{The estimates p-value.}
#' \item{bootstrapped}{Bootstrapped auxiliary statistics.}}
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom zeallot %<-%
#'
#' @export
KPSS.N.breaks.bootstrap <- function(y, x,
                                    model, break.point,
                                    const = FALSE, trend = FALSE,
                                    weakly.exog = TRUE,
                                    lags.init, leads.init,
                                    corr.max, kernel,
                                    iter = 9999,
                                    bootstrap = "sample",
                                    criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
    }

    N <- nrow(y)

    c(., test, u, ., DOLS.lags, DOLS.leads, .) %<-%
        KPSS.N.breaks(
            y, x,
            model, break.point,
            const, trend,
            weakly.exog,
            lags.init, leads.init,
            corr.max, kernel,
            criterion
        )

    if (weakly.exog) {
        xreg <- cbind(
            x,
            determinants.KPSS.N.breaks(model, N, break.point, const, trend)
        )
    } else {
        c(., xreg) %<-%
            DOLS.vars.N.breaks(
                y, x,
                model, break.point,
                const, trend,
                DOLS.lags, DOLS.leads
            )
    }

    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    registerDoSNOW(cluster)

    result <- foreach(
        i = 1:iter,
        .combine = rbind,
        .options.snow = list(progress = progress)
    ) %dopar% {
        if (bootstrap == "sample") {
            temp.y <- sample(u, length(u), replace = TRUE)
        } else if (bootstrap == "Cavaliere-Taylor") {
            z <- rnorm(length(u))
            temp.y <- z * u
        } else if (bootstrap == "Rademacher") {
            z <- sample(c(-1, 1), length(u), replace = TRUE)
            temp.y <- z * u
        }

        c(beta, resid, ., t.beta) %<-% OLS(temp.y, xreg)

        S.t <- apply(resid, 2, cumsum)
        if (!is.null(kernel)) {
            temp.test <- N^(-2) * drop(t(S.t) %*% S.t) /
                alrvr.kernel(resid, corr.max, kernel)
        } else {
            temp.test <- N^(-2) * drop(t(S.t) %*% S.t) / alrvr(resid)
        }
        temp.test
    }

    stopCluster(cluster)

    p.value <- (1 / iter) * sum(I(test <= result))

    return(
        list(
            test = test,
            p.value = p.value,
            bootstrapped = result
        )
    )
}
