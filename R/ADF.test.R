#' @title
#' ADF.test - simple implementation of ADF test
#'
#' @param y A variable of interest.
#' @param const Include const to the model if TRUE.
#' @param trend Include trend to the model if TRUE.
#' @param max.lag Maximum lag number
#' @param criterion A criterion used to select number of lags.
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
                     criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (! criterion %in% c("bic", "aic", "lwz")) {
        warning("WARNING! Unknown criterion, BIC is used")
        criterion <- "bic"
    }

    N <- nrow(y)
    pos <- ifelse(const && trend, 3, ifelse(const, 2, 1))

    d.y <- diff(y)
    x <- y[1:(N - 1), , drop = FALSE]

    deter <- NULL
    if (const)
        deter <- cbind(deter, rep(1, (N - 1)))
    if (trend)
        deter <- cbind(deter, 1:(N - 1))

    if (max.lag > 0) {
        for (l in 1:max.lag) {
            x <- cbind(x, lagn(d.y, l))
        }
    }

    c(min.beta, min.resid, ., min.t.beta) %<-%
        OLS(d.y, cbind(deter, x[, 1]))
    min.lag <- 0

    if (max.lag > 0) {
        if (criterion == "bic")
            min.criterion <- log(drop(t(min.resid) %*% min.resid) /
                                 (length(min.resid) - length(min.beta))) +
                length(min.beta) * log(length(min.resid) - length(min.beta)) /
                (length(min.resid) - length(min.beta))
        else if (criterion == "aic")
            min.criterion <- log(drop(t(min.resid) %*% min.resid) /
                                 (length(min.resid) - length(min.beta))) +
                2 * length(min.beta) / (length(min.resid) - length(min.beta))
        else if (criterion == "lwz")
            min.criterion <- log(drop(t(min.resid) %*% min.resid) /
                                 (length(min.resid) - length(min.beta))) +
                0.299 * length(min.beta) *
                (log(length(min.resid) - length(min.beta)))^2.1

        for (l in 1:max.lag) {
            deter <- NULL
            if (const)
                deter <- cbind(deter, rep(1, (N - 1 - l)))
            if (trend)
                deter <- cbind(deter, 1:(N - 1 - l))
            c(tmp.beta, tmp.resid, ., tmp.t.beta) %<-%
                OLS(
                    d.y[(l + 1):(N - 1), , drop = FALSE],
                    cbind(
                        deter,
                        x[(l + 1):(N - 1), 1:(l + 1), drop = FALSE]
                    )
                )
            if (criterion == "bic")
                tmp.criterion <- log(drop(t(tmp.resid) %*% tmp.resid) /
                                     (length(tmp.resid) - length(tmp.beta))) +
                    length(tmp.beta) * log(length(tmp.resid) - length(tmp.beta)) /
                    (length(tmp.resid) - length(tmp.beta))
            else if (criterion == "aic")
                tmp.criterion <- log(drop(t(tmp.resid) %*% tmp.resid) /
                                     (length(tmp.resid) - length(tmp.beta))) +
                    2 * length(tmp.beta) / (length(tmp.resid) - length(tmp.beta))
            else if (criterion == "lwz")
                tmp.criterion <- log(drop(t(tmp.resid) %*% tmp.resid) /
                                     (length(tmp.resid) - length(tmp.beta))) +
                    0.299 * length(tmp.beta) *
                    (log(length(tmp.resid) - length(tmp.beta)))^2.1

            if (tmp.criterion < min.criterion) {
                min.criterion <- tmp.criterion
                min.beta <- tmp.beta
                min.resid <- tmp.resid
                min.t.beta <- tmp.t.beta
                min.lag <- l
            }
        }
    }

    return(
        list(
            y = drop(y),
            const = const,
            trend = trend,
            beta = min.beta,
            t.beta = drop(min.t.beta[pos]),
            critical.value = critical.value,
            residuals = min.resid,
            lag = min.lag
        )
    )
}
