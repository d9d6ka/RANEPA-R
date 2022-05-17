#' @title
#' Produce lagged vector backward of forward
#'
#' @param x Initial vector.
#' @param i Size of lag (lead if negative).
#'
#' @return Lagged or leaded vector.
lagn <- function(x, i) {
    if (!is.matrix(x)) x <- as.matrix(x)
    N <- nrow(x)
    k <- ncol(x)
    if (i > 0)
        return(
            rbind(
                matrix(data = NA, nrow = i, ncol = k),
                x[1:(N - i), , drop = FALSE]
            )
        )
    else
        return(
            rbind(
                x[(1 + abs(i)):N, , drop = FALSE],
                matrix(data = NA, nrow = abs(i), ncol = k)
            )
        )
}
