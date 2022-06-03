seasonal.dummies <- function(N) {
    s1 <- c(1 - 1 / 12, rep(-1 / 12, 11))

    result <- NULL
    for (i in 0:10) {
        result <- cbind(
            result,
            c(
                rep(-1 / 12, i),
                rep(s1, length.out = N - i)
            )
        )
    }

    return(result)
}
