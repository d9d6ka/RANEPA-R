olsqr <- function(y, x) {
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    p <- x %*% b
    r <- y - p
    return(
        list(beta = b, resid = r, predict = p)
    )
}
