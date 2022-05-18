#' @importFrom zeallot %<-%
ols.residuals.N.breaks <- function(y, x, break.point) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y)
    ntb <- length(break.point)

    break.point <- c(0, break.point, N)

    result <- NULL

    for (i in 1:(ntb + 1)) {
        temp.y <- y[(break.point[i] + 1):break.point[i + 1], , drop = FALSE]
        temp.x <- x[(break.point[i] + 1):break.point[i + 1], , drop = FALSE]
        c(., r, ., .) %<-% OLS(temp.y, temp.x)
        if (i == 1) result <- r
        else result <- rbind(result, r)
    }

    return(result)
}
