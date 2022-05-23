#' @title
#' NW.estimation - Nadaraya–Watson kernel regression.
#'
#' @details
#' See Harvey, Leybourne and Zu (2021) for more information.
#'
#' @param y LHS dependent variable.
#' @param x RHS explanation variable.
#' @param h Bandwidth.
#' @param kernel Needed kernel, currently only `unif` and `gauss`.
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
        if (kernel == "unif") {
            W <- ifelse(abs((k - (1:N)) / N / h) <= 1, 1, 0)
        } else if (kernel == "gauss") {
            W <- pnorm((k - (1:N)) / N / h)
        }
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
