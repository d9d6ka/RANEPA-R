#' @title
#' Custom OLS with extra information
#'
#' @description
#' Getting OLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `resid`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`.
#'
#' @importFrom stats .lm.fit
#'
#' @keywords internal
OLS <- function(y,
                x) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    tmp.model <- .lm.fit(x, y)
    S.2 <- drop(t(tmp.model$residuals) %*% tmp.model$residuals) /
        (nrow(x) - ncol(x))
    t.beta <- tmp.model$coefficients / sqrt(diag(S.2 * qr.solve(t(x) %*% x)))

    return(
        list(
            beta = as.matrix(tmp.model$coefficients),
            residuals = tmp.model$residuals,
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
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param c A coefficient for \eqn{\rho} calculation.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `resid`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`.
#'
#' @keywords internal
GLS <- function(y,
                z,
                c) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(z)) z <- as.matrix(z)

    n.obs <- nrow(y)
    n.var <- ncol(y)

    rho <- 1 + c / n.obs

    y.hat <- y - rho * lagn(y, 1)
    y.hat[1, ] <- y[1, ]

    z.hat <- z - rho * lagn(z, 1)
    z.hat[1, ] <- z[1, ]

    betas <- NULL
    resids <- NULL
    fitted.values <- NULL
    t.betas <- NULL

    for (i in 1:n.var) {
        res.OLS <- OLS(y.hat[, i, drop = FALSE], z.hat)
        betas <- cbind(betas, res.OLS$beta)
        t.betas <- cbind(
            t.betas,
            res.OLS$t.beta
        )
        fitted.values <- cbind(
            fitted.values,
            z %*% res.OLS$beta
        )
        resids <- cbind(
            resids,
            y[, i, drop = FALSE] - fitted.values
        )
    }

    return(
        list(
            beta = betas,
            residuals = resids,
            predict = fitted.values,
            t.beta = drop(t.betas)
        )
    )
}


#' @title
#' Custom AR with extra information
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param max.lag The maximum number of lags.
#' @param criterion A criterion for lag number estimation.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `residuals`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `lag`: estimated number of lags.
#'
#' @keywords internal
AR <- function(y,
               x,
               max.lag,
               criterion = "aic") {
    if (!is.null(criterion)) {
        if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
            warning("WARNING! Unknown criterion, none is used")
            criterion <- NULL
        }
    }

    if (!is.matrix(y)) y <- as.matrix(y)
    n.obs <- nrow(y)
    tmp.y <- y[(1 + max.lag):n.obs, , drop = FALSE]

    if (!is.null(x)) {
        if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)
        k <- ncol(x)
        tmp.x <- x[(1 + max.lag):n.obs, , drop = FALSE]
    } else {
        k <- 0
        tmp.x <- NULL
    }

    for (l in 1:max.lag) {
        if (l <= max.lag) {
            tmp.x <- cbind(
                tmp.x,
                lagn(y, l)[(1 + max.lag):n.obs, , drop = FALSE]
            )
        }
    }

    if (is.null(criterion)) {
        res.lag <- max.lag
        tmp.OLS <- OLS(tmp.y, tmp.x[, 1:(k + res.lag), drop = FALSE])
        res.beta <- tmp.OLS$beta
        res.resid <- tmp.OLS$residuals
        res.predict <- tmp.OLS$predict
        res.t.beta <- tmp.OLS$t.beta
        rm(tmp.OLS)
    } else {
        res.lag <- 0

        if (!is.null(x)) {
            tmp.OLS <- OLS(tmp.y, tmp.x[, 1:k, drop = FALSE])
            res.beta <- tmp.OLS$beta
            res.resid <- tmp.OLS$residuals
            res.predict <- tmp.OLS$predict
            res.t.beta <- tmp.OLS$t.beta
            rm(tmp.OLS)
            res.IC <- log(drop(t(res.resid) %*% res.resid) / (n.obs - max.lag))
        } else {
            res.IC <- Inf
        }

        for (l in 1:max.lag) {
            if (l <= max.lag) {
                tmp.OLS <- OLS(tmp.y, tmp.x[, 1:(k + l), drop = FALSE])
                tmp.beta <- tmp.OLS$beta
                tmp.resid <- tmp.OLS$residuals
                tmp.predict <- tmp.OLS$predict
                tmp.t.beta <- tmp.OLS$t.beta
                rm(tmp.OLS)
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
    }

    return(
        list(
            beta = res.beta,
            residuals = res.resid,
            predict = res.predict,
            t.beta = res.t.beta,
            lag = res.lag,
            criterion = res.IC
        )
    )
}


#' @title
#' Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated coefficient vector and
#' residuals.
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.estimation <- function(y,
                          x,
                          h,
                          kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    n.obs <- length(y)

    rho <- rep(0, n.obs)
    for (k in 1:n.obs) {
        W <- NW.kernel(k, (1:n.obs) / n.obs, h, kernel)
        rho[k] <- sum(x * W * y) / sum(x * W * x)
    }

    return(
        list(
            my = y,
            mx = x,
            h = h,
            kernel = kernel,
            rr1.est = rho,
            u.hat = y - rho * x
        )
    )
}


#' @title
#' Nadaraya–Watson kernel volatility estimation
#'
#' @param e A series of interest.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated omega and s.e.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor.
#' “Lag Length Selection for Unit Root Tests in the Presence
#' of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.volatility <- function(e,
                          h,
                          kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    n.obs <- length(e)

    omega.sq <- rep(0, n.obs)
    for (k in 1:n.obs) {
        W <- NW.kernel(k, (1:n.obs) / n.obs, h, kernel)
        omega.sq[k] <- sum(W * e^2) / sum(W)
    }

    return(
        list(
            me = e,
            h = h,
            kernel = kernel,
            omega.sq = omega.sq,
            se = sqrt(omega.sq)
        )
    )
}


#' @title
#' LOO-CV for h in Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x An explanatory variable.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @return A list of arguments as well as the estimated bandwidth `h`.
#'
#' @keywords internal
NW.loocv <- function(y,
                     x,
                     kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    n.obs <- length(y)

    HT <- seq(n.obs^(-0.5), n.obs^(-0.3), by = 0.01)
    cv0 <- Inf

    for (hi in HT) {
        rho <- rep(0, n.obs)
        for (k in 1:n.obs) {
            W <- NW.kernel(k, (1:n.obs) / n.obs, hi, kernel)
            W[k] <- 0
            rho[k] <- sum(x * W * y) / sum(x * W * x)
        }

        cv1 <- sum((y - rho * x)^2)
        if (cv1 < cv0) {
            cv0 <- cv1
            h <- hi
        }
    }

    return(
        list(
            my = y,
            mx = x,
            kernel = kernel,
            h = h
        )
    )
}


#' @title
#' Returning a kernel value for a specific point
#'
#' @param i An index of the base point.
#' @param x A series for kernel calculations.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
NW.kernel <- function(i,
                      x,
                      h,
                      kernel = "unif") {
    if (kernel == "unif") {
        W <- ifelse(abs((x - x[i]) / h) <= 1, 1, 0)
    } else if (kernel == "gauss") {
        W <- pnorm((x - x[i]) / h)
    }
    return(W)
}
