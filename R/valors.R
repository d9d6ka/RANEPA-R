#' @export
critical_values_kpss_1p <- function(model, break.point, N, k) { # nolint
    lambda <- break.point / N
    i <- trunc(lambda * 10)
    if (i == 0) i <- 1
    if (i == 10) i <- 9

    m_vc <- .cval_kpss_1p[[model]][[k]]

    return(m_vc[, i, drop = FALSE])
}

#' @export
critical_values_kpss_2p <- function(model, break.point, N) {
    lambda1 <- round(break.point[1] / N, 1)
    i <- trunc(lambda1 * 10)

    lambda2 <- round(break.point[2] / N, 1)
    j <- trunc(lambda2 * 10)
    if (1 <= model & model <= 4)
        j <- j - 1

    v_vc <- .cval_kpss_2p[[model]]

    return(
        c(
            v_vc[1][i, j],
            v_vc[2][i, j],
            v_vc[3][i, j],
            v_vc[4][i, j]
        )
    )
}
