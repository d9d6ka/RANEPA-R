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
ssr_partition_mp <- function(y, x, m = 1, width = 2, ssr_data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    if (is.null(ssr_data))
        ssr_data <- ssr_matrix(y, x, width)

    if (m == 1) {
        tmp_result <- ssr_partition_1p(
            1, N,
            width, N - width + 1,
            N, ssr_data
        )
        optimal_ssr <- tmp_result$ssr
        optimal_break <- tmp_result$break_point
    }
    else {
        variants <- N - (m + 1) * width + 1
        loop_ssr <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        loop_break <- matrix(
            data = 0,
            nrow = variants,
            ncol = m
        )
        for (step in 1:m) {
            tmp_ssr <- matrix(
                data = Inf,
                nrow = variants,
                ncol = 1
            )
            if (step == 1) {
                for (v in 1:variants) {
                    loop_end <- 2 * width + v - 1
                    tmp_res <- ssr_partition_1p(
                        1,
                        loop_end,
                        width,
                        loop_end - width,
                        N,
                        ssr_data
                    )
                    loop_ssr[v, 1] <- tmp_res$ssr
                    loop_break[v, 1] <- tmp_res$break_point
                }
            }
            else if (step == m) {
                for (v in 1:variants) {
                    tmp_ssr[v, 1] <- loop_ssr[v, 1] +
                        ssr_data[step * width + v, N]
                }
                optimal_ssr <- min(tmp_ssr)
                optimal_index <- which.min(tmp_ssr)
                optimal_break <- loop_break[optimal_index, ]
                optimal_break[m] <- step * width + optimal_index - 1
            }
            else {
                next_ssr <- matrix(
                    data = Inf,
                    nrow = variants,
                    ncol = 1
                )
                next_break <- matrix(
                    data = 0,
                    nrow = variants,
                    ncol = m
                )
                for (loop_end in ((step + 1) * width):(N - (m - step) * width)) { # nolint
                    next_v <- loop_end - (step + 1) * width + 1
                    for (v in 1:variants) {
                        tmp_ssr[v, 1] <- loop_ssr[v, 1] +
                            ssr_data[step * width + v, loop_end]
                    }
                    next_ssr[next_v, 1] <- min(tmp_ssr)
                    next_index <- which.min(tmp_ssr)
                    next_break[next_v, 1:m] <- loop_break[next_index, ]
                    next_break[next_v, step] <- step * width + next_index - 1
                }
                loop_ssr <- next_ssr
                loop_break <- next_break
            }
        }
    }

    return(
        list(
            ssr = optimal_ssr,
            break_point = optimal_break
        )
    )
}
