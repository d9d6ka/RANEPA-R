determi_kpss_1p <- function(model, t, tb) {
    du <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1, nrow = t - tb, ncol = 1)
    )
    dt <- rbind(
        matrix(data = 0, nrow = tb, ncol = 1),
        matrix(data = 1:(t - tb), nrow = t - tb, ncol = 1)
    )

    if (model == 1) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du
        )
    }
    else if (model == 2) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1)
        )
    }
    else if (model == 3) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else if (model == 4) {
        xt <- cbind(
            matrix(data = 1, nrow = t, ncol = 1),
            du,
            matrix(data = 1:t, nrow = t, ncol = 1),
            dt
        )
    }
    else {
        stop("ERROR: Try to specify the deterministic component again")
    }

    return(xt)
}
