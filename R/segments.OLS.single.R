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
#' \item{break_point}{The point of possible break.}
#' }
segments.OLS.single <- function(beg, end,
                                first.break, last.break,
                                len, SSR.data) {
    tmp.result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (tb in first.break:last.break) {
        tmp.result[tb] <- SSR.data[beg, tb] +
            SSR.data[tb + 1, end]
    }

    min.ssr <- min(tmp.result[first.break:last.break])
    break.point <- (first.break - 1) +
        which.min(tmp.result[first.break:last.break])

    return(
        list(
            SSR = min.ssr,
            break.point = break.point
        )
    )
}
