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
    tmp.model <- lm.fit(x, y)

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
