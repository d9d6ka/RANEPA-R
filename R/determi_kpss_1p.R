#' @title
#' Construct determinant variables.
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 1 structural break.
#'
#' @param model \describe{
#' \item{1}{Model with trend, break in const.}
#' \item{2}{Model with const and trend, break in const.}
#' \item{3}{Model with const and trend, break in trend.}
#' \item{4}{Model with const and trend, break in const and trend.}
#' }
#' @param N Number of observations.
#' @param break.point Break point.
#'
#' @return Matrix of determinant variables.
#'
determi_kpss_1p <- function(model, N, break.point) { # nolint
    du <- rbind(
        matrix(data = 0, nrow = break.point, ncol = 1),
        matrix(data = 1, nrow = N - break.point, ncol = 1)
    )
    dt <- rbind(
        matrix(data = 0, nrow = break.point, ncol = 1),
        matrix(data = 1:(N - break.point), nrow = N - break.point, ncol = 1)
    )
    const <- matrix(data = 1, nrow = N, ncol = 1)
    trend <- matrix(data = 1:N, nrow = N, ncol = 1)

    if (model == 1) {
        xt <- cbind(
            trend,
            du
        )
    }
    else if (model == 2) {
        xt <- cbind(
            const,
            du,
            trend
        )
    }
    else if (model == 3) {
        xt <- cbind(
            const,
            trend,
            dt
        )
    }
    else if (model == 4) {
        xt <- cbind(
            const,
            du,
            trend,
            dt
        )
    }
    else {
        stop("ERROR: Try to specify the deterministic component again")
    }

    return(xt)
}
