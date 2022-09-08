#' @title
#' Estimating optimal number of leads and lags
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y LHS dependent variable.
#' @param trend Whether the trend is to be included.
#' @param zb I(1) regressors with break.
#' @param zf I(1) regressors without break.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break date.
#' @param criterion A criterion for lead and lag number estimation.
#'
#' @return A list of estimated values of leads and lags.
#'
#' @references
#' Kurozumi, Eiji, and Anton Skrobotov.
#' “Confidence Sets for the Break Date in Cointegrating Regressions.”
#' Oxford Bulletin of Economics and Statistics 80, no. 3 (2018): 514–35.
#' https://doi.org/10.1111/obes.12223.
select.lead.lag <- function(y,
                            trend = TRUE,
                            zb = NULL,
                            zf = NULL,
                            trim = 0.05,
                            criterion = "bic") {
    if (!criterion %in% c("bic", "aic", "hq", "lwz")) {
        stop("ERROR! Unknown criterion")
    }

    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    first.break <- trunc(2 * trim * n.obs)
    last.break <- trunc((1 - 2 * trim) * n.obs)

    wb <- cbind(
        rep(1, n.obs),
        if (trend) (1:n.obs) / n.obs else NULL,
        zb
    )

    z <- cbind(zb, zf)
    w <- cbind(wb, zf)

    u.hat <- OLS(y, w)$residuals
    ssr.0 <- drop(t(u.hat) %*% u.hat)
    est.date <- n.obs

    for (t in first.break:last.break) {
        wb1 <- rbind(
            matrix(data = 0, nrow = t, ncol = 1),
            as.matrix(wb[(t + 1):n.obs, ])
        )
        w <- cbind(wb, wb1, zf)
        u.hat <- OLS(y, w)$residuals
        ssr.1 <- drop(t(u.hat) %*% u.hat)

        if (ssr.1 < ssr.0) {
            ssr.0 <- ssr.1
            est.date <- t
        }
    }

    est.dt <- rbind(
        matrix(data = 0, nrow = est.date, ncol = ncol(wb)),
        as.matrix(wb[(est.date + 1):n.obs, ])
    )
    y <- as.matrix(y[2:n.obs, ])
    wb <- as.matrix(wb[2:n.obs, ])
    est.dt <- as.matrix(est.dt[2:n.obs, ])
    d.z <- as.matrix(z[2:n.obs, ]) - as.matrix(z[1:(n.obs - 1), ])

    wf <- cbind(
        if (!is.null(zf)) as.matrix(zf[2:n.obs, ]) else NULL,
        d.z
    )

    N <- nrow(y)

    max.lead.lag <- trunc(4 * (N / 100)^(1 / 4))
    if ((est.date - max.lead.lag - 1) <= ncol(wb)) {
        max.lead.lag <- est.date - ncol(wb) - 2
    } else if ((N - max.lead.lag - est.date) <= ncol(wb)) {
        max.lead.lag <- N - est.date - ncol(wb) - 1
    }

    y.0 <- as.matrix(y[(max.lead.lag + 1):(N - max.lead.lag), ])
    w.0 <- cbind(
        as.matrix(wb[(max.lead.lag + 1):(N - max.lead.lag), ]),
        as.matrix(est.dt[(max.lead.lag + 1):(N - max.lead.lag), ]),
        as.matrix(wf[(max.lead.lag + 1):(N - max.lead.lag), ])
    )

    u.hat <- OLS(y.0, w.0)$residuals
    min.ic <- info.criterion(u.hat, ncol(w.0))[[criterion]]
    est.lead <- 0
    est.lag <- 0

    for (cur.lead in 1:max.lead.lag) {
        for (cur.lag in 1:max.lead.lag) {
            w.1 <- w.0
            for (k in 1:cur.lead) {
                w.1 <- cbind(
                    w.1,
                    as.matrix(
                        d.z[(max.lead.lag + 1 - k):(N - max.lead.lag - k), ]
                    )
                )
            }
            for (k in 1:cur.lag) {
                w.1 <- cbind(
                    w.1,
                    as.matrix(
                        d.z[(max.lead.lag + 1 + k):(N - max.lead.lag + k), ]
                    )
                )
            }

            u.hat <- OLS(y.0, w.1)$residuals

            cur.ic <- info.criterion(u.hat, ncol(w.1))[[criterion]]
            if (cur.ic < min.ic) {
                est.lead <- cur.lead
                est.lag <- cur.lag
                min.ic <- cur.ic
            }
        }
    }

    return(
        list(
            lead = est.lead,
            lag = est.lag
        )
    )
}
