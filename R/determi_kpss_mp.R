## model = 1 - const
## model = 2 - trend
## model = 3 - const + trend
determi_kpss_mp <- function(model, N, break_point, trend = FALSE) { # nolint
    ntb <- length(break_point)
    if (length(model) == 1)
        model <- rep(model, ntb)
    else if (length(model) != ntb)
        stop("ERROR! Inconsistent sizes of `model` and `break_point`")

    xt <- matrix(data = 1, nrow = N, ncol = 1)
    if (trend)
        xt <- cbind(
            xt,
            matrix(data = 1:N, nrow = N, ncol = 1)
        )

    for (i in 1:ntb)
        if (model[i] == 1)
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break_point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = N - break_point[i],
                        ncol = 1
                    )
                )
            )
        else if (model[i] == 2)
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break_point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(N - break_point[i]),
                        nrow = N - break_point[i],
                        ncol = 1
                    )
                )
            )
        else if (model[i] == 3)
            xt <- cbind(
                xt,
                rbind(
                    matrix(
                        data = 0,
                        nrow = break_point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1,
                        nrow = N - break_point[i],
                        ncol = 1
                    )
                ),
                rbind(
                    matrix(
                        data = 0,
                        nrow = break_point[i],
                        ncol = 1
                    ),
                    matrix(
                        data = 1:(N - break_point[i]),
                        nrow = N - break_point[i],
                        ncol = 1
                    )
                )
            )

    return(xt)
}
