lagn <- function(x, i) {
    if (!is.matrix(x)) x <- as.matrix(x)
    t <- nrow(x)
    k <- ncol(x)
    if (i > 0)
        return(
            rbind(
                matrix(data = NA, nrow = i, ncol = k),
                x[1:(t - i), , drop = FALSE]
            )
        )
    else
        return(
            rbind(
                x[(1 + abs(i)):t, , drop = FALSE],
                matrix(data = NA, nrow = abs(i), ncol = k)
            )
        )
}
