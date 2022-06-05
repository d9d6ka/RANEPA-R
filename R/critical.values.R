critical.values.KPSS.1p <- function(model, break.point, N, k) {
    lambda <- round(break.point / N, 1)
    i <- trunc(lambda * 10)
    if (i == 0) i <- 1
    if (i == 10) i <- 9

    m_vc <- .cval_KPSS_1p[[model]][[k]]

    return(m_vc[, i, drop = FALSE])
}

critical.values.KPSS.2p <- function(model, break.point, N) {
    lambda1 <- round(break.point[1] / N, 1)
    i <- trunc(lambda1 * 10)

    lambda2 <- round(break.point[2] / N, 1)
    j <- trunc(lambda2 * 10)
    if (1 <= model & model <= 4) {
        j <- j - 1
    }

    v_vc <- .cval_KPSS_2p[[model]]

    return(
        c(
            v_vc[1][i, j],
            v_vc[2][i, j],
            v_vc[3][i, j],
            v_vc[4][i, j]
        )
    )
}

p.values.SADF <- function(statistic, N.obs, cr.values) {
    N.list <- as.numeric(names(cr.values))

    if (N.obs < min(N.list)) {
        warning("Too little number of observations, using data for T = 30")
        i.0 <- min(N.list)
    } else {
        i.0 <- max(N.list[N.list <= N.obs])
    }
    p.0 <- sum(cr.values[[as.character(i.0)]] > statistic) /
        length(cr.values[[as.character(i.0)]])

    if (N.obs > max(N.list)) {
        i.1 <- max(N.list)
    } else {
        i.1 <- min(N.list[N.list >= N.obs])
    }
    p.1 <- sum(cr.values[[as.character(i.1)]] > statistic) /
        length(cr.values[[as.character(i.1)]])

    if (i.0 != i.1) {
        p.value <- p.0 + (p.1 - p.0) * (N - i.0) / (i.1 - i.0)
    } else {
        p.value <- p.0
    }
    return(p.value)
}
