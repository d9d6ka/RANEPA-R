#' Procedure to compute deterministic terms
#' for KPSS with 2 structural breaks.
#'
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#' @param N Number of observations.
#' @param tb1 Position for the first structural break
#'            (respective to the origin which is 1).
#' @param tb2 Position for the second structural break
#'            (respective to the origin which is 1).
#'
#' @return Matrix of deterministic terms.
determi_kpss_2p <- function(model, N, tb1, tb2) {
    cns <- matrix(data = 1, nrow = N, ncol = 1)
    tend <- matrix(data = 1:N, nrow = N, ncol = 1)

    if (model == 1) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, du1, du2)
    }
    else if (model == 2) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, du2)
    }
    else if (model == 3) {
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1:(N - tb2), nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, dt1, dt2)
    }
    else if (model == 4) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        dt2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1:(N - tb2), nrow = N - tb2, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2, dt2)
    }
    else if (model == 5) {
        du1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1, nrow = N - tb1, ncol = 1)
        )
        du2 <- rbind(
            matrix(data = 0, nrow = tb2, ncol = 1),
            matrix(data = 1, nrow = N - tb2, ncol = 1)
        )
        dt1 <- rbind(
            matrix(data = 0, nrow = tb1, ncol = 1),
            matrix(data = 1:(N - tb1), nrow = N - tb1, ncol = 1)
        )
        xt <- cbind(cns, tend, du1, dt1, du2)
    }
    else
        stop("Try to speciy another model")

    return(xt)
}
