#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A dependent (LHS) variable.
#' @param x A matrix of explanatory (RHS) variables.
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param break.point An array of moments of structural breaks.
#' @param const,trend Whether a constant or trend are to be included.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of
#' \itemize{
#' \item Estimates of coefficients,
#' \item Estimates of residuals,
#' \item A set of informational criterions values,
#' \item \eqn{t}-statistics for the estimates of coefficients.}
#'
#' @importFrom zeallot %<-%
DOLS.N.breaks <- function(y,
                          x,
                          model,
                          break.point,
                          const = FALSE,
                          trend = FALSE,
                          k.lags,
                          k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) {
        stop("ERROR! Explanatory variables needed for DOLS")
    }
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    c(yreg, xreg) %<-%
        DOLS.vars.N.breaks(
            y, x,
            model, break.point,
            const, trend,
            k.lags, k.leads
        )

    c(beta, resid, ., t.beta) %<-% OLS(yreg, xreg)

    s2 <- drop(t(resid) %*% resid) / nrow(xreg)

    criterions <- info.criterion(resid, ncol(xreg))

    return(
        list(
            beta = beta,
            resid = resid,
            criterions = criterions,
            t.beta = t.beta
        )
    )
}


#' @title
#' Preparing variables for DOLS regression with multiple known break points
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A dependent (LHS) variable.
#' @param x A matrix of explanatory (RHS) variables.
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param break.point An array of moments of structural breaks.
#' @param const,trend Whether a constant or trend are to be included.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of LHS and RHS variables.
DOLS.vars.N.breaks <- function(y, x,
                               model, break.point,
                               const = FALSE, trend = FALSE,
                               k.lags, k.leads) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(x)) stop("ERROR! Explanatory variables needed for DOLS")
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    d.x.step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE]
    d.x.lag <- d.x.step
    d.x.lead <- d.x.step

    for (i in 1:k.lags) {
        d_x <- cbind(
            d.x.lag,
            lagn(d.x.step, i)
        )
    }

    for (i in 1:k.leads) {
        d.x.lead <- cbind(
            d.x.lead,
            lagn(d.x.step, -i)
        )
    }

    if (k.lags != 0 & k.leads != 0) {
        lags <- d.x.lag
        leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
        lags.leads <- cbind(lags, leads)
        lags.leads <- lags.leads[(k.lags + 1):(N - 1 - k.leads), , drop = FALSE]
    } else if (k.lags != 0 & k.leads == 0) {
        lags <- d.x.lag
        lags.leads <- lags[(k.lags + 1):(N - 1), , drop = FALSE]
    } else if (k.lags == 0 & k.leads != 0) {
        lags <- d.x.lag
        leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
        lags.leads <- cbind(lags, leads)
        lags.leads <- lags.leads[1:(N - 1 - k.leads), , drop = FALSE]
    } else if (k.lags == 0 & k.leads == 0) {
        lags.leads <- d.x.lag
    }
    deter <- determinants.KPSS.N.breaks(model, N, break.point, const, trend)

    xreg <- cbind(
        deter[(k.lags + 2):(N - k.leads), , drop = FALSE],
        x[(k.lags + 2):(N - k.leads), , drop = FALSE],
        lags.leads
    )

    yreg <- y[(k.lags + 2):(N - k.leads), 1, drop = FALSE]

    return(
        list(
            yreg = yreg,
            xreg = xreg
        )
    )
}
