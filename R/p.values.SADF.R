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