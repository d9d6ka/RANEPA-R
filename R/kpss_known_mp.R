#' @importFrom zeallot %<-%
#' 
#' @export
kpss_known_mp <- function(y, model, break_point, kmax, kernel, trend = FALSE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y) # nolint

    z <- determi_kpss_mp(model, N, break_point, trend)
    c(beta, resid, p, t_b) %<-% olsqr(y, z)

    s_t <- apply(resid, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr_kernel(resid, kmax, kernel)
    else
        test <- N^(-2) * drop(t(s_t) %*% s_t) / alrvr(resid)

    return(
        list(
            test = test,
            break_point = break_point
        )
    )
}
