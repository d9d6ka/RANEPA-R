#' @title
#' Construct determinant variables.
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 2 structural breaks.
#'
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' \item{6}{for the AC-CA model.}
#' \item{7}{for the AC-CA model.}
#' }
#' @param N Number of observations.
#' @param break.point Positions for the first and second structural breaks
#'            (respective to the origin which is 1).
#'
#' @return Matrix of deterministic terms.
determinants.KPSS.2.breaks <- function(model, N, break.point) {
    cns <- matrix(data = 1, nrow = N, ncol = 1)
    tend <- matrix(data = 1:N, nrow = N, ncol = 1)

    if (model == 1) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, du1, du2)
    }
    else if (model == 2) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, du2)
    }
    else if (model == 3) {
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1:(N - break.point[1]),
                   nrow = N - break.point[1], ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1:(N - break.point[2]),
                   nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, dt1, dt2)
    }
    else if (model == 4) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = N - break.point[2], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1:(N - break.point[1]),
                   nrow = N - break.point[1], ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1:(N - break.point[2]),
                   nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2, dt2)
    }
    else if (model == 5) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1:(N - break.point[2]),
                   nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt2)
    }
    else if (model == 6) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1, nrow = N - break.point[2], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1:(N - break.point[1]),
                   nrow = N - break.point[1], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2)
    }
    else if (model == 7) {
        du1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1, nrow = N - break.point[1], ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = break.point[1], ncol = 1),
            matrix(data = 1:(N - break.point[1]),
                   nrow = N - break.point[1], ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = break.point[2], ncol = 1),
            matrix(data = 1:(N - break.point[2]),
                   nrow = N - break.point[2], ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, dt2)
    }
    else
        stop("Try to speciy another model")

    return(xt)
}
