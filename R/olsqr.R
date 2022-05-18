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
olsqr <- function(y, x) {
    beta <- qr.solve(t(x) %*% x) %*% t(x) %*% y
    predict <- x %*% beta
    resid <- y - predict
    s2 <- drop(t(resid) %*% resid) / (nrow(x) - ncol(x))
    t.beta <- sweep(beta, 1, sqrt(diag(s2 * qr.solve(t(x) %*% x))), `/`)
    return(
        list(
            beta = beta,
            resid = resid,
            predict = predict,
            t.beta = t.beta
        )
    )
}
