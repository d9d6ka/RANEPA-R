#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @details
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param beg Sample begin.
#' @param end Sample end.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param len Total number of observations.
#' @param SSR.data The matrix of recursive SSR values.
#'
#' @return List containing
#' \describe{
#' \item{SSR}{Optimal SSR value.}
#' \item{break_point}{The breaking point.}
#' }
segs.SSR.1.break <- function(beg, end, first.break, last.break, len, SSR.data) {
    tmp_result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (p_break in first.break:last.break) {
        tmp_result[p_break] <- SSR.data[beg, p_break] +
            SSR.data[p_break + 1, end]
    }

    min.SSR <- min(tmp_result[first.break:last.break])
    break.point <- (first.break - 1) +
        which.min(tmp_result[first.break:last.break])

    return(
        list(
            SSR = min.SSR,
            break.point = break.point
        )
    )
}
