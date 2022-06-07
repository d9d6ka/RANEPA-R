#' @title
#' ADF.test - simple implementation of ADF test
#'
#' @param y A variable of interest.
#' @param const Include const to the model if TRUE.
#' @param trend Include trend to the model if TRUE.
#' @param max.lag Maximum lag number
#' @param criterion A criterion used to select number of lags.
#' If lag selection is not needed keep this NULL.
#' @param modified.criterion Whether the unit-root test modificaton is needed.
#'
#' @return List containing:
#' \itemize{
#' \item y
#' \item const
#' \item trend
#' \item residuals
#' \item coefficient estimates
#' \item t-statistic value
#' \item critical value
#' \item Number of lags
#' \item indicator of stationarity}
#'
#' @importFrom zeallot %<-%
#' @export
ADF.test <- function(y, const = TRUE, trend = FALSE, max.lag = 0,
                     criterion = NULL, modified.criterion = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.null(criterion)) {
        if (criterion %in% c("bic", "aic", "lwz", "hq")) {
            warning("WARNING! Unknown criterion, none is used")
            criterion <- NULL
        }
    }

    N <- nrow(y)
    pos <- ifelse(const && trend, 3, ifelse(const, 2, 1))

    deter <- NULL
    if (const) {
        deter <- cbind(deter, rep(1, N - 1))
    }
    if (trend) {
        deter <- cbind(deter, 1:(N - 1))
    }

    d.y <- diff(y)
    x <- cbind(
        deter,
        y[1:(N - 1), , drop = FALSE]
    )

    if (max.lag > 0) {
        for (l in 1:max.lag) {
            x <- cbind(x, lagn(d.y, l))
        }
    }

    if (is.null(criterion)) {
        res.lag <- max.lag
    } else {
        c(., e, ., .) %<-% OLS(
            d.y[(1 + max.lag):(N - 1), , drop = FALSE],
            x[(1 + max.lag):(N - 1), 1:pos, drop = FALSE]
        )

        res.IC <- log(drop(t(e) %*% e) / N)
        res.lag <- 0

        for (l in 1:max.lag) {
            if (max.lag == 0) break
            c(., e, ., .) %<-% OLS(
                d.y[(1 + max.lag):(N - 1), , drop = FALSE],
                x[(1 + max.lag):(N - 1), 1:(pos + l), drop = FALSE]
            )
            tmp.IC <- info.criterion(e, l)[[criterion]]
            if (tmp.IC < res.IC) {
                res.IC <- tmp.IC
                res.lag <- l
            }
        }
    }

    c(res.beta, res.resid, ., res.t.beta) %<-%
        OLS(
            d.y[(max.lag + 1):(N - 1), , drop = FALSE],
            x[(max.lag + 1):(N - 1), 1:(pos + res.lag), drop = FALSE]
        )

    return(
        list(
            y = drop(y),
            const = const,
            trend = trend,
            beta = res.beta,
            t.beta = drop(res.t.beta),
            alpha = res.beta[pos],
            t.alpha = drop(res.t.beta[pos]),
            residuals = res.resid,
            lag = res.lag
        )
    )
}
