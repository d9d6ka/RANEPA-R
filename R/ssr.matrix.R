#' @title
#' Pre-calculate matrix of recursive SSR values.
#'
#' @param y Dependent variable.
#' @param x Explanatory variables.
#' @param width Minimum spacing between the breaks.
#'
#' @return The matrix of recursive SSR values.
#'
#' @export
ssr.matrix <- function(y, x, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    result <- matrix(data = Inf, nrow = N, ncol = N)

    for (i in 1:(N - width + 1)) {
        result[i, 1:N] <- ssr.recursive(
            y, x, i, N, width
        )
    }

    return(result)
}
