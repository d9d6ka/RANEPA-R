kpss_unknown_mp <- function(y, x, m = 1, width = 2, ssr_data = NULL) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    N <- nrow(y) # nolint

    if (is.null(ssr_data))
        ssr_data <- ssr_matrix(y, x, width)
}
