#' @title
#' Deterministic terms for multi-break KPSS
#'
#' @description
#' Procedure to compute deterministic terms
#' for KPSS with \eqn{m} structural breaks.
#'
#' @details
#' **model** should be either a scalar or a vector of the same size as the
#' **break.point**. If scalar **model** will be repeated till the length of
#' **break.point** is achieved.
#'
#' @param model A scalar or vector of
#' \describe{
#' \item{1}{for the break in const.}
#' \item{2}{for the break in trend.}
#' \item{3}{for the break in const and trend.}
#' }
#' @param N Number of observations.
#' @param break.point Array of structural breaks.
#' @param const Include constant if **TRUE**.
#' @param trend Include trend if **TRUE**.
#'
#' @return Matrix of deterministic terms.
determinants.KPSS.N.breaks <- function(model, N, break.point,
                                       const = FALSE, trend = FALSE) {
    ntb <- length(break.point)
    if (length(model) == 1) {
        model <- rep(model, ntb)
    } else if (length(model) != ntb) {
        stop("ERROR! Inconsistent sizes of model and break.point")
    }

    xt <- NULL
    if (const) {
        xt <- matrix(data = 1, nrow = N, ncol = 1)
    }
    if (trend) {
        xt <- cbind(
            xt,
            matrix(data = 1:N, nrow = N, ncol = 1)
        )
    }

    for (i in 1:ntb) {
        if (model[i] == 1) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = N - break.point[i],
                        ncol = 1
                    )
                )
            )
        } else if (model[i] == 2) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(N - break.point[i]),
                        nrow = N - break.point[i],
                        ncol = 1
                    )
                )
            )
        } else if (model[i] == 3) {
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = N - break.point[i],
                        ncol = 1
                    )
                ),
                rbind(
                    matrix(
                        data = 0,
                        nrow = break.point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(N - break.point[i]),
                        nrow = N - break.point[i],
                        ncol = 1
                    )
                )
            )
        }
    }

    return(xt)
}
