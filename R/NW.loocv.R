#' @title
#' NW.loocv - LOO-CV for h in Nadaraya–Watson kernel regression.
#'
#' @details
#' See Harvey, Leybourne and Zu (2021) for more information.
#'
#' @param y LHS dependent variable.
#' @param x RHS explanation variable.
#' @param kernel Needed kernel, currently only `unif` and `gauss`.
#'
#' @return A list of arguments as well as the estimated bandwidth `h`.
NW.loocv <- function(y, x, kernel = "unif") {
    if (!kernel %in% c("unif", "gauss")) {
        warning("WARNING! Unknown kernel, unif is used instead")
        kernel <- "unif"
    }

    N <- length(y)

    HT <- seq(N ^ (-0.5), N ^ (-0.3), by = 0.01)
    cv0 <- Inf

    for (hi in seq_len(HT)) {
        h1 <- HT[hi]

        rho <- rep(0, N)
        for (k in 1:N) {
            if (kernel == "unif") {
                W <- ifelse(abs(((1:N) - k) / N / h1) <= 1, 1, 0)
            } else if (kernel == "gauss") {
                W <- pnorm(((1:N) - k) / N / h1)
            }
            W[k] <- 0
            rho[k] <- sum(x * W * y) / sum(x * W * x)
        }

        cv1 <- sum((y - rho * x)^2)
        if (cv1 < cv0) {
            cv0 <- cv1
            h <- h1
        }
    }

    return(
        list(
            my = y,
            mx = x,
            kernel = kernel,
            h.est = h
        )
    )
}
