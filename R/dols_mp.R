#' @import MASS
dols_mp <- function(y, x,
                    model, break.point, # nolint
                    const = FALSE, trend = FALSE,
                    k.lags, k.leads) { # nolint
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    c(yreg, xreg) %<-%
        dols_prepare_mp(y, x, model, break.point, const, trend, k.lags, k.leads)

    c(beta, resid, ., t_beta) %<-% olsqr(yreg, xreg)

    s2 <- drop(t(resid) %*% resid) / (nrow(xreg) - ncol(xreg))

    bic <- log(s2) + ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

    return(
        list(
            beta   = beta,
            resid  = resid,
            bic    = bic,
            t_beta = t_beta
        )
    )
}
