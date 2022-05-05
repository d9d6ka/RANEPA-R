ssr_matrix <- function(y, x, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)

    result <- matrix(data = Inf, nrow = N, ncol = N)

    for (i in 1:(N - width + 1)) {
        result[i, 1:N] <- ssr_recursive(
            y, x, i, N, width
        )
    }

    return(result)
}
