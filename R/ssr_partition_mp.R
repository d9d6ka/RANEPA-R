#' Find m + 1 optimal pertitions
#'
#' @details
#' Based on Bai and Perron (2003).
#'
#' @details
#' Based on Brown, Durbin and Evans (1975).
#' 
#' @param y (Tx1)-vector of the dependent variable
#' @param x (Txk)-vector of the explanatory stochastic regressors
#' @param m Number of breaks
#' @param width Minimum spacing between the breaks
#' 
#' @return List of 2 elements: optimal SSR and the vector of breaks
#' 
#' @export
ssr_partition_mp <- function(y, x, m = 1, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)
    ssr_data <- ssr_matrix(y, x, width)

    if (m == 1) {
        tmp_result <- ssr_partition_1p(
            1, t,
            width, t - width + 1,
            t, ssr_data
        )
        optimal_ssr <- tmp_result$ssr
        optimal_break <- tmp_result$break_point
    }
    else {
        variants <- t - (m + 1) * width + 1
        loop_ssr <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        loop_ends <- matrix(
            data = 0,
            nrow = variants,
            ncol = 1
        )
        loop_break <- matrix(
            data = 0,
            nrow = variants,
            ncol = m
        )
        tmp_ssr <- matrix(
            data = Inf,
            nrow = variants,
            ncol = 1
        )
        for (step in 1:m) {
            if (step == 1) {
                for (v in 1:variants) {
                    loop_end <- 2 * width + v - 1
                    tmp_res <- ssr_partition_1p(
                        1,
                        loop_end,
                        width,
                        loop_end - width,
                        t,
                        ssr_data
                    )
                    loop_ssr[v, 1] <- tmp_res$ssr
                    loop_break[v, 1] <- tmp_res$break_point
                    loop_ends[v, 1] <- loop_end
                }
            }
            else if (step == m) {
                for (v in 1:variants) {
                    tmp_ssr[v, 1] <- loop_ssr[v, 1] +
                        ssr_data[step * width + v, t]
                }
                optimal_ssr <- min(tmp_ssr)
                optimal_index <- which.min(tmp_ssr)
                optimal_break <- loop_break[optimal_index, ]
                optimal_break[m] <- step * width + optimal_index
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
                for (loop_end in ((step + 1) * width):(t - (m - step) * width)) { # nolint
                    next_v <- loop_end - (step + 1) * width + 1
                    for (v in 1:variants) {
                        tmp_ssr[v, 1] <- loop_ssr[v, 1] +
                            ssr_data[step * width + v, loop_end]
                    }
                    next_ssr[next_v, 1] <- min(tmp_ssr)
                    next_index <- which.min(tmp_ssr)
                    next_break[next_v, 1:m] <- loop_break[next_index, ]
                    next_break[next_v, step] <- step * width + next_index
                    loop_ends[next_v, 1] <- loop_end
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
