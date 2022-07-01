#' @title
#' KPSS-test with multiple unknown structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with multiple unknown structural breaks
#'
#' @param y An input (LHS) time series of interest.
#' @param x A matrix of (RHS) explanatory stochastic regressors.
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param break.point Array of structural breaks.
#' @param const Include constant if **TRUE**.
#' @param trend Include trend if **TRUE**.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param lags.init,leads.init Scalars defininig the initial number of lags and
#' leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the BIC
#' information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#' @param iter Number of bootstrap iterations.
#' @param bootstrap Type of bootstrapping:
#' * `"sample"`: sampling from residuals with replacement,
#' * `"Cavaliere-Taylor"`: multiplying residuals by \eqn{N(0, 1)}-distributed
#' variable,
#' * `"Rademacher"`: multiplying residuals by Rademacher-distributed variable.
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic or lwz.
#'
#' @return A list of:
#' * `test`: The value of KPSS test statistic,
#' * `p.value`: The estimates p-value,
#' * `bootstrapped`: Bootstrapped auxiliary statistics.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
KPSS.N.breaks.bootstrap <- function(y, x,
                                    model, break.point,
                                    const = FALSE, trend = FALSE,
                                    weakly.exog = TRUE,
                                    lags.init, leads.init,
                                    max.lag, kernel,
                                    iter = 9999,
                                    bootstrap = "sample",
                                    criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
    }

    N <- nrow(y)

    tmp.kpss <- KPSS.N.breaks(
        y, x,
        model, break.point,
        const, trend,
        weakly.exog,
        lags.init, leads.init,
        max.lag, kernel,
        criterion
    )
    test <- tmp.kpss$test
    u <- tmp.kpss$residuals
    DOLS.lags <- tmp.kpss$DOLS.lags
    DOLS.leads <- tmp.kpss$DOLS.leads
    rm(tmp.kpss)

    if (weakly.exog) {
        xreg <- cbind(
            x,
            determinants.KPSS.N.breaks(model, N, break.point, const, trend)
        )
    } else {
        xreg <- DOLS.vars.N.breaks(
            y, x,
            model, break.point,
            const, trend,
            DOLS.lags, DOLS.leads
        )$xreg
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

        resids <- OLS(temp.y, xreg)$residuals

        if (!is.null(kernel)) {
            temp.test <- KPSS(resids, lr.var.SPC(resids, max.lag, kernel))
        } else {
            temp.test <- KPSS(resids, lr.var.bartlett.AK(resids))
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
