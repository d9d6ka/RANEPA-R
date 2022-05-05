#' @importFrom zeallot %<-%
resid_mp <- function(y, x, break_point) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    break_point <- c(0, break_point)

    result <- NULL

    for (i in seq_len(break_point)) {
        tmp_y <- y[(break_point[i] + 1):break_point[i + 1], , drop = FALSE]
        tmp_x <- x[(break_point[i] + 1):break_point[i + 1], , drop = FALSE]
        c(beta, r, p, t_b) %<-% olsqr(tmp_y, tmp_x)
        if (i == 1) result <- r
        else result <- rbind(result, r)
    }

    return(result)
}
