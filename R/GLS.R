#' @importFrom zeallot %<-%
GLS <- function(y, z, c) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(z)) z <- as.matrix(z)

    N <- nrow(y)

    rho <- 1 + c / N

    y.hat <- y - rho * lagn(y, 1)
    y.hat[1, ] <- y[1, ]

    z.hat <- z - rho * lagn(z, 1)
    z.hat[1, ] <- z[1, ]

    c(beta, ., ., t.beta) %<-% OLS(y.hat, z.hat)
    predict <- y - z %*% beta
    resid <- y - predict

    return(
        list(
            beta = beta,
            resid = resid,
            predict = predict,
            t.beta = t.beta
        )
    )
}
