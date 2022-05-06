#' Procedure to minimize the SSR for 2 break points
#'
#' @param beg Sample begin.
#' @param end Sample end.
#' @param first_b First possible break point.
#' @param last_b Last possible break point.
#' @param len Total number of observations.
#' @param ssr_data The matrix of recursive SSR values.
#'
#' @return List containing
#' \describe{
#' \item{ssr}{Optimal SSR value.}
#' \item{break_point}{The breaking point.}
#' }
ssr_partition_1p <- function(beg, end, first_b, last_b, len, ssr_data) {
    tmp_result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (p_break in first_b:last_b) {
        tmp_result[p_break] <- ssr_data[beg, p_break] +
            ssr_data[p_break + 1, end]
    }

    tmp_ssr <- min(tmp_result[first_b:last_b])
    tmp_break <- (first_b - 1) + which.min(tmp_result[first_b:last_b])

    return(
        list(
            ssr = tmp_ssr,
            break_point = tmp_break
        )
    )
}
