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

    N <- nrow(y)

    first.break <- trunc(2 * trim * N)
    last.break <- trunc((1 - 2 * trim) * N)

    wb <- cbind(
        rep(1, N),
        if (trend) (1:N) / N else NULL,
        zb
    )

    z <- cbind(zb, zf)

    w <- cbind(wb, zf)

    b.hat <- solve(t(w) %*% t) %*% t(w) %*% y
    u.hat <- y - w %*% b.hat
    ssr.0 <- drop(t(u.hat) %*% u.hat)
    est.date <- N

    for (t in first.break:last.break) {
        wb1 <- rbind(
            matrix(data = 0, nrow = t, ncol = 1),
            as.matrix(wb[(t + 1):N, ])
        )
        w <- cbind(wb, wb1, zf)
        b.hat <- solve(t(w) %*% t) %*% t(w) %*% y
        u.hat <- y - w %*% b.hat
        ssr.1 <- drop(t(u.hat) %*% u.hat)

        if (ssr.1 < ssr.0) {
            ssr.0 <- ssr.1
            est.date <- t
        }
    }

    est.dt <- rbind(
        matrix(data = 0, nrow = est.date, ncol = ncol(wb)),
        as.matrix(wb[(est.date + 1):N, ])
    )
    y <- as.matrix(y[2:N, ])
    wb <- as.matrix(wb[2:N, ])
    est.dt <- as.matrix(est.dt[2:N, ])
    d.z <- as.matrix(z[2:N, ]) - as.matrix(z[1:(N - 1), ])

    wf <- cbind(
        if (!is.null(zf)) as.matrix(zf[2:N, ]) else NULL,
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

    b.hat <- solve(t(w.0) %*% w.0) %*% t(w.0) %*% y.0
    u.hat <- y.0 - w.0 %*% b.hat
    min.ic <- info.criterion(u.hat, ncol(w.0))[[criterion]]
    est.lead <- 0
    est.lag <- 0

    for (cur.lead in 1:max.lead.lag) {
        for (cur.lag in 1:max.lead.lag) {
            w_1 <- w.0
            for (k in 1:cur.lead) {
                w_1 <- cbind(
                    w_1,
                    as.matrix(
                        d.z[(max.lead.lag + 1 - k):(N - max.lead.lag - k), ]
                    )
                )
            }
            for (k in 1:cur.lag) {
                w_1 <- cbind(
                    w_1,
                    as.matrix(
                        d.z[(max.lead.lag + 1 + k):(N - max.lead.lag + k), ]
                    )
                )
            }

            b.hat <- solve(t(w_1) %*% w_1) %*% t(w_1) %*% y.0
            u.hat <- y.0 - w_1 %*% b.hat

            cur.ic <- info.criterion(u.hat, ncol(w_1))[[criterion]]
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
