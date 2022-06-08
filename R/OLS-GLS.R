#' @title
#' Custom OLS with extra information
#'
#' @description
#' Getting OLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y Dependent variable.
#' @param x Explanatory variables.
#'
#' @return The list of betas, residuals, forecasted values and t-values.
#'
#' @import MASS
OLS <- function(y, x) {
    tmp.model <- .lm.fit(x, y)
    S.2 <- drop(t(tmp.model$residuals) %*% tmp.model$residuals) /
        (nrow(x) - ncol(x))
    t.beta <- tmp.model$coefficients / sqrt(diag(S.2 * qr.solve(t(x) %*% x)))

    return(
        list(
            beta = tmp.model$coefficients,
            resid = tmp.model$residuals,
            predict = tmp.model$fitted.values,
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
#' @param y Dependent variable.
#' @param x Explanatory variables.
#' @param c Coefficient for \eqn{\rho} calculation.
#'
#' @return The list of betas, residuals, forecasted values and t-values.
#'
#' @import MASS
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
