#' @importFrom zeallot %<-%
kpss_unknown_mp <- function(y, x, m = 1, width = 2, ssr_data = NULL, kmax = 0, kernel = "bartlett") {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    if (is.null(ssr_data))
        ssr_data <- ssr_matrix(y, x, width)

    c(ssr, tb) %<-% ssr_partition_mp(y, x, m, width, ssr_data)
    resid <- resid_mp(y, x, tb)

    s_t <- apply(resid, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr_kernel(resid, kmax, kernel)
    else
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)

    return(
        list(
            test = test,
            break_point = tb
        )
    )
}
