ssr_matrix <- function(y, x, width = 2) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    t <- nrow(y)

    result <- matrix(data = Inf, nrow = t, ncol = t)

    for (i in 1:(t - width + 1)) {
        result[i, 1:t] <- ssr_recursive(
            y, x, i, t, width
        )
    }

    return(result)
}
