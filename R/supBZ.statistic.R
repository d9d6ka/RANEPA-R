#' @title
#' Calculate supBZ statistic
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y The series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param sigma.sq Local non-parametric estimates of variance. If `NULL` they
#' will be estimated via Nadaraya-Watson procedure.
#' @param generalized Whether to calculate generalized statistic value.
#'
#' @return A list of:
#' * `y`,
#' * `trim`,
#' * `sigma.sq`,
#' * `BZ.values`: a series of BZ-statistic,
#' * `supBZ.value`: the maximum of `supBZ.values`,
#' * `h.est`: the estimated value of bandwidth if `sigma.sq` is `NULL`.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Testing Explosive Bubbles with Time-Varying Volatility.”
#' Econometric Reviews 38, no. 10 (November 26, 2019): 1131–51.
#' https://doi.org/10.1080/07474938.2018.1536099.
#'
#' @export
supBZ.statistic <- function(y,
                            trim = 0.01 + 1.8 / sqrt(length(y)),
                            sigma.sq = NULL,
                            generalized = FALSE) {
    n.obs <- length(y)

    if (is.null(sigma.sq)) {
        ## NW estimation.
        my <- (diff(y))^2
        mx <- rep(1, n.obs - 1)
        nw.loocv.model <- NW.loocv(my, mx, kernel = "gauss")
        h.est <- nw.loocv.model$h
        nw.model <- NW.volatility(
            my,
            kernel = "gauss",
            h = nw.loocv.model$h
        )
        sigma.sq <- nw.model$omega.sq
    }

    y <- y - y[1]
    d.y <- diff(y)
    l.y <- y[1:(n.obs - 1)]

    BZ.values <- c()
    m <- 1

    if (!generalized) {
        for (j in (floor(trim * n.obs)):n.obs) {
            BZ.values[m] <-
                sum(d.y[1:(j - 1)] * l.y[1:(j - 1)] / sigma.sq[1:(j - 1)]) /
                    (sum(l.y[1:(j - 1)]^2 / sigma.sq[1:(j - 1)]))^0.5
            m <- m + 1
        }
    } else {
        for (i in 1:(n.obs - floor(trim * n.obs) + 1)) {
            for (j in (i + floor(trim * n.obs) - 1):n.obs) {
                BZ.values[m] <-
                    sum(d.y[i:(j - 1)] * l.y[i:(j - 1)] / sigma.sq[i:(j - 1)]) /
                        (sum(l.y[i:(j - 1)]^2 / sigma.sq[i:(j - 1)]))^0.5
                m <- m + 1
            }
        }
    }

    supBZ.value <- max(BZ.values)

    return(
        c(
            list(
                y = y,
                trim = trim,
                sigma.sq = sigma.sq,
                BZ.values = BZ.values,
                supBZ.value = supBZ.value
            ),
            if (exists("h.est")) {
                list(h.est = h.est)
            } else {
                NULL
            }
        )
    )
}
