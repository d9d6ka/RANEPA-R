#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param beg Sample begin.
#' @param end Sample end.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param len Total number of observations.
#' @param SSR.data The matrix of recursive SSR values.
#'
#' @return A list of:
#' * `SSR`: Optimal SSR value,
#' * `break.point`: The point of possible break.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
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
