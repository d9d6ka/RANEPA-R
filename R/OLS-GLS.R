#' @title
#' Custom OLS with extra information
#'
#' @description
#' Getting OLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y Dependent variable.
#' @param x Explanatory variables.
#'
#' @return The list of betas, residuals, forecasted values and t-values.
OLS <- function(y, x) {
    tmp.model <- .lm.fit(x, y)
    fitted.values <- y - tmp.model$residuals
    S.2 <- drop(t(tmp.model$residuals) %*% tmp.model$residuals) /
        (nrow(x) - ncol(x))
    t.beta <- tmp.model$coefficients / sqrt(diag(S.2 * qr.solve(t(x) %*% x)))

    return(
        list(
            beta = as.matrix(tmp.model$coefficients),
            resid = tmp.model$residuals,
            predict = fitted.values,
            t.beta = t.beta
        )
    )
}


#' @title
#' Custom GLS with extra information
#'
#' @description
#' Getting GLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y Dependent variable.
#' @param z Explanatory variables.
#' @param c Coefficient for \eqn{\rho} calculation.
#'
#' @return The list of betas, residuals, forecasted values and t-values.
#'
#' @importFrom zeallot %<-%
GLS <- function(y, z, c) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(z)) z <- as.matrix(z)

    N <- nrow(y)

    rho <- 1 + c / N

    y.hat <- y - rho * lagn(y, 1)
    y.hat[1, ] <- y[1, ]

    z.hat <- z - rho * lagn(z, 1)
    z.hat[1, ] <- z[1, ]

    c(beta, ., ., t.beta) %<-% OLS(y.hat, z.hat)
    predict <- z %*% beta
    resid <- y - predict

    return(
        list(
            beta = beta,
            resid = resid,
            predict = predict,
            t.beta = t.beta
        )
    )
}


#' @title
#' Custom AR with extra information
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y Dependent variable.
#' @param x Exogenous explanatory variables.
#' @param max.lag The maximum number of lags.
#' @param criterion A criterion for lag number estimation.
#'
#' @importFrom zeallot %<-%
AR <- function(y, x, max.lag, criterion = "aic") {
    if (!is.null(criterion)) {
        if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
            warning("WARNING! Unknown criterion, none is used")
            criterion <- NULL
        }
    }

    if (!is.matrix(y)) y <- as.matrix(y)
    N <- nrow(y)
    tmp.y <- y[(1 + max.lag):N, , drop = FALSE]

    if (!is.null(x)) {
        if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)
        k <- ncol(x)
        tmp.x <- x[(1 + max.lag):N, , drop = FALSE]
    } else {
        k <- 0
        tmp.x <- NULL
    }

    for (l in 1:max.lag) {
        tmp.x <- cbind(
            tmp.x,
            lagn(y, l)[(1 + max.lag):N, , drop = FALSE]
        )
    }

    if (is.null(criterion)) {
        res.lag <- max.lag
        c(res.beta, res.resid, res.predict, res.t.beta) %<-%
            OLS(tmp.y, tmp.x[, 1:(k + res.lag), drop = FALSE])
    } else {
        res.lag <- 0

        if (!is.null(x)) {
            c(res.beta, res.resid, res.predict, res.t.beta) %<-%
                OLS(tmp.y, tmp.x[, 1:k, drop = FALSE])
            res.IC <- log(drop(t(resid) %*% resid) / (N - max.lag))
        } else {
            res.IC <- Inf
        }

        for (l in 1:max.lag) {
            c(tmp.beta, tmp.resid, tmp.predict, tmp.t.beta) %<-%
                OLS(tmp.y, tmp.x[, 1:(k + l), drop = FALSE])
            temp.IC <- info.criterion(tmp.resid, l)[[criterion]]

            if (temp.IC < res.IC) {
                res.IC <- temp.IC
                res.beta <- tmp.beta
                res.resid <- tmp.resid
                res.predict <- tmp.predict
                res.t.beta <- tmp.t.beta
                res.lag <- l
            }
        }
    }

    return(
        list(
            beta = res.beta,
            resid = res.resid,
            predict = res.predict,
            t.beta = res.t.beta,
            lag = res.lag
        )
    )
}
