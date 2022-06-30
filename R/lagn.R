#' @title
#' Produce a vector lagged backward of forward
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param x Initial vector.
#' @param i Size of lag (lead if negative).
#' @param na Value to fill missing observations, `NA` by default.
#'
#' @return Lagged or leaded vector.
lagn <- function(x, i, na = NA) {
    if (!is.matrix(x)) x <- as.matrix(x)
    N <- nrow(x)
    k <- ncol(x)
    if (i > 0) {
        return(
            rbind(
                matrix(data = na, nrow = i, ncol = k),
                x[1:(N - i), , drop = FALSE]
            )
        )
    } else {
        return(
            rbind(
                x[(1 + abs(i)):N, , drop = FALSE],
                matrix(data = na, nrow = abs(i), ncol = k)
            )
        )
    }
}
