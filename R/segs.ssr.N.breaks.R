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
segs.ssr.N.breaks <- function(y, x, m = 1, width = 2, ssr.data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    if (is.null(ssr.data))
        ssr.data <- ssr.matrix(y, x, width)

    if (m == 1) {
        tmp.result <- segs.ssr.1.break(
            1, N,
            width, N - width + 1,
            N, ssr.data
        )
        optimal.ssr <- tmp.result$ssr
        optimal.break <- tmp.result$break.point
    }
    else {
        variants <- N - (m + 1) * width + 1
        step.ssr <- matrix(
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
            temp.ssr <- matrix(
                data = Inf,
                nrow = variants,
                ncol = 1
            )
            if (step == 1) {
                for (v in 1:variants) {
                    step.end <- 2 * width + v - 1
                    tmp_res <- segs.ssr.1.break(
                        1,
                        step.end,
                        width,
                        step.end - width,
                        N,
                        ssr.data
                    )
                    step.ssr[v, 1] <- tmp_res$ssr
                    step.break[v, 1] <- tmp_res$break.point
                }
            }
            else if (step == m) {
                for (v in 1:variants) {
                    temp.ssr[v, 1] <- step.ssr[v, 1] +
                        ssr.data[step * width + v, N]
                }
                optimal.ssr <- min(temp.ssr)
                optimal.index <- which.min(temp.ssr)
                optimal.break <- step.break[optimal.index, ]
                optimal.break[m] <- step * width + optimal.index - 1
            }
            else {
                next.ssr <- matrix(
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
                        temp.ssr[v, 1] <- step.ssr[v, 1] +
                            ssr.data[step * width + v, step.end]
                    }
                    next.ssr[next.v, 1] <- min(temp.ssr)
                    next.index <- which.min(temp.ssr)
                    next.break[next.v, 1:m] <- step.break[next.index, ]
                    next.break[next.v, step] <- step * width + next.index - 1
                }
                step.ssr <- next.ssr
                step.break <- next.break
            }
        }
    }

    return(
        list(
            ssr = optimal.ssr,
            break.point = optimal.break
        )
    )
}
