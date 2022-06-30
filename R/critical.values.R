#' @title
#' Auxiliary function returning pre-calculated critical values for
#' [KPSS.1.break()].
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param model \describe{
#' \item{1}{for model An.}
#' \item{2}{for model A.}
#' \item{3}{for model B.}
#' \item{4}{for model C.}
#' \item{5}{for model D.}
#' \item{6}{for model E.}
#' }
#' @param break.point Position of the break point.
#' @param N Number of observations.
#' @param k Number of RHS variables.
critical.values.KPSS.1p <- function(model, break.point, N, k) {
    lambda <- round(break.point / N, 1)
    i <- trunc(lambda * 10)
    if (i == 0) i <- 1
    if (i == 10) i <- 9

    m_vc <- .cval_KPSS_1p[[model]][[k]]

    return(m_vc[, i, drop = FALSE])
}

#' @title
#' Auxiliary function returning pre-calculated critical values for
#' [KPSS.2.breaks()].
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#' @param break.point Position of the break point.
#' @param N Number of observations.
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


#' @title
#' Interpolating p-value for intermediate observation numbers for SADF-type
#' tests.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param statistic The statistic value.
#' @param N.obs The number of observations.
#' @param cr.values The set of precalculated tables.
p.values.SADF <- function(statistic, N.obs, cr.values) {
    N.table.obs <- as.numeric(names(cr.values))

    if (N.obs < min(N.table.obs)) {
        warning("Too little number of observations, using data for T = 30")
        i.0 <- min(N.table.obs)
    } else {
        i.0 <- max(N.table.obs[N.table.obs <= N.obs])
    }
    p.0 <- sum(cr.values[[as.character(i.0)]] > statistic) /
        length(cr.values[[as.character(i.0)]])

    if (N.obs > max(N.table.obs)) {
        i.1 <- max(N.table.obs)
    } else {
        i.1 <- min(N.table.obs[N.table.obs >= N.obs])
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
