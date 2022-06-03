#' @export
info.criterion <- function(resid, extra, criterion, ...) {
    if (!is.matrix(resid)) resid <- as.matrix(resid)

    N <- nrow(resid)

    if (criterion == "bic") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            extra * log(N - extra) / (N - extra)
    } else if (criterion == "aic") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            2 * extra / (N - extra)
    } else if (criterion == "lwz") {
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            0.299 * extra * (log(N - extra))^2.1
    } else if (criterion == "maic") {
        s2 <- drop(t(resid) %*% resid) / (N - extra)
        tau <- (...elt(1) ^ 2) * drop(t(...elt(2)) %*% ...elt(2)) / s2
        result <- log(drop(t(resid) %*% resid) / (N - extra)) +
            2 * (tau + extra) / (N - extra)
    }

    return(result)
}
