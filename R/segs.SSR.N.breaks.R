#' @title
#' Find \eqn{m + 1} optimal partitions
#'
#' @details
#' Based on Bai and Perron (2003).
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#'
#' @return List of 2 elements: optimal SSR and the vector of breaks.
#'
#' @export
segs.SSR.N.breaks <- function(y, x, m = 1, width = 2, SSR.data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    if (is.null(SSR.data)) {
        SSR.data <- SSR.matrix(y, x, width)
    }

    if (m == 1) {
        tmp.result <- segs.SSR.1.break(
            1, N,
            width, N - width + 1,
            N, SSR.data
        )
        optimal.SSR <- tmp.result$SSR
        optimal.break <- tmp.result$break.point
    } else {
        variants <- N - (m + 1) * width + 1
        step.SSR <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        step.break <- matrix(
            data = 0,
            nrow = variants,
            ncol = m
        )
        for (step in 1:m) {
            temp.SSR <- matrix(
                data = Inf,
                nrow = variants,
                ncol = 1
            )
            if (step == 1) {
                for (v in 1:variants) {
                    step.end <- 2 * width + v - 1
                    tmp_res <- segs.SSR.1.break(
                        1,
                        step.end,
                        width,
                        step.end - width,
                        N,
                        SSR.data
                    )
                    step.SSR[v, 1] <- tmp_res$SSR
                    step.break[v, 1] <- tmp_res$break.point
                }
            } else if (step == m) {
                for (v in 1:variants) {
                    temp.SSR[v, 1] <- step.SSR[v, 1] +
                        SSR.data[step * width + v, N]
                }
                optimal.SSR <- min(temp.SSR)
                optimal.index <- which.min(temp.SSR)
                optimal.break <- step.break[optimal.index, ]
                optimal.break[m] <- step * width + optimal.index - 1
            } else {
                next.SSR <- matrix(
                    data = Inf,
                    nrow = variants,
                    ncol = 1
                )
                next.break <- matrix(
                    data = 0,
                    nrow = variants,
                    ncol = m
                )
                for (step.end in ((step + 1) * width):(N - (m - step) * width)) {
                    next.v <- step.end - (step + 1) * width + 1
                    for (v in 1:variants) {
                        temp.SSR[v, 1] <- step.SSR[v, 1] +
                            SSR.data[step * width + v, step.end]
                    }
                    next.SSR[next.v, 1] <- min(temp.SSR)
                    next.index <- which.min(temp.SSR)
                    next.break[next.v, 1:m] <- step.break[next.index, ]
                    next.break[next.v, step] <- step * width + next.index - 1
                }
                step.SSR <- next.SSR
                step.break <- next.break
            }
        }
    }

    return(
        list(
            SSR = optimal.SSR,
            break.point = optimal.break
        )
    )
}
