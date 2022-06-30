#' @title
#' NW.estimation - Nadaraya–Watson kernel regression.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y LHS dependent variable.
#' @param x RHS explanation variable.
#' @param h Bandwidth.
#' @param kernel Needed kernel, currently only `unif` and `gauss`.
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @return A list of arguments as well as the estimated coefficient vector and
#' residuals.
NW.estimation <- function(y, x, h, kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    N <- length(y)

    rho <- rep(0, N)
    for (k in 1:N) {
        W <- NW.kernel(k, (1:N) / N, h, kernel)
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
#' NW.volatility - Nadaraya–Watson kernel volatility estimation
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param e The series of interest.
#' @param h Bandwidth.
#' @param kernel Needed kernel, currently only `unif` and `gauss`.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor.
#' “Lag Length Selection for Unit Root Tests in the Presence
#' of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' @return A list of arguments as well as the estimated omega and s.e.
NW.volatility <- function(e, h, kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    N <- length(e)

    omega.sq <- rep(0, N)
    for (k in 1:N) {
        W <- NW.kernel(k, (1:N) / N, h, kernel)
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
#' NW.loocv - LOO-CV for h in Nadaraya–Watson kernel regression.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y LHS dependent variable.
#' @param x RHS explanation variable.
#' @param kernel Needed kernel, currently only `unif` and `gauss`.
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @return A list of arguments as well as the estimated bandwidth `h`.
NW.loocv <- function(y, x, kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    N <- length(y)

    HT <- seq(N^(-0.5), N^(-0.3), by = 0.01)
    cv0 <- Inf

    for (hi in HT) {
        rho <- rep(0, N)
        for (k in 1:N) {
            W <- NW.kernel(k, (1:N) / N, hi, kernel)
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

NW.kernel <- function(i, x, h, kernel = "unif") {
    if (kernel == "unif") {
        W <- ifelse(abs((x - x[i]) / h) <= 1, 1, 0)
    } else if (kernel == "gauss") {
        W <- pnorm((x - x[i]) / h)
    }
    return(W)
}
