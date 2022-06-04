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
    if (!(is.null(criterion) || criterion %in% c("bic", "aic", "lwz"))) {
        warning("WARNING! Unknown criterion, none is used")
        criterion <- NULL
    }

    aux.deter <- function(const, trend, lag, obs) {
        deter <- NULL
        if (const) {
            deter <- cbind(deter, rep(1, (N - 1 - lag)))
        }
        if (trend) {
            deter <- cbind(deter, 1:(N - 1 - lag))
        }
        return(deter)
    }

    N <- nrow(y)
    pos <- ifelse(const && trend, 3, ifelse(const, 2, 1))

    d.y <- diff(y)
    x <- y[1:(N - 1), , drop = FALSE]

    if (max.lag > 0) {
        for (l in 1:max.lag) {
            x <- cbind(x, lagn(d.y, l))
        }
    }

    c(min.beta, min.resid, ., min.t.beta) %<-%
        OLS(
            d.y[(max.lag + 1):(N - 1), , drop = FALSE],
            cbind(
                aux.deter(const, trend, max.lag, N),
                x[(max.lag + 1):(N - 1), 1:(max.lag + 1), drop = FALSE]
            )
        )
    min.lag <- max.lag

    if (max.lag > 0 && !is.null(criterion)) {
        min.criterion <- info.criterion(
            min.resid,
            length(min.beta),
            criterion,
            modified.criterion,
            min.beta[pos],
            x[(max.lag + 1):(N - 1), 1, drop = FALSE]
        )

        for (l in seq(max.lag - 1, 0, -1)) {
            deter <- NULL
            if (const) {
                deter <- cbind(deter, rep(1, (N - 1 - l)))
            }
            if (trend) {
                deter <- cbind(deter, 1:(N - 1 - l))
            }
            c(tmp.beta, tmp.resid, ., tmp.t.beta) %<-%
                OLS(
                    d.y[(l + 1):(N - 1), , drop = FALSE],
                    cbind(
                        deter,
                        x[(l + 1):(N - 1), 1:(l + 1), drop = FALSE]
                    )
                )
            tmp.criterion <- info.criterion(
                tmp.resid,
                length(tmp.beta),
                criterion,
                modified.criterion,
                tmp.beta[pos],
                x[(l + 1):(N - 1), 1, drop = FALSE]
            )

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
            residuals = min.resid,
            lag = min.lag
        )
    )
}
