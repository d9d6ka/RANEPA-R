#' @import MASS
DOLS.N.breaks <- function(y, x,
                    model, break.point,
                    const = FALSE, trend = FALSE,
                    k.lags, k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    c(yreg, xreg) %<-%
        DOLS.vars.N.breaks(y, x,
                           model, break.point,
                           const, trend,
                           k.lags, k.leads)

    c(beta, resid, ., t.beta) %<-% OLS(yreg, xreg)

    s2 <- drop(t(resid) %*% resid) / (nrow(xreg) - ncol(xreg))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta   = beta,
            resid  = resid,
            bic    = bic,
            t.beta = t.beta
        )
    )
}
