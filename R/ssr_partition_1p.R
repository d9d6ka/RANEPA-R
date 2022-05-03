ssr_partition_1p <- function(beg, end, first_b, last_b, len, ssr_data) {
    tmp_result <- matrix(data = Inf, nrow = len, ncol = 1)

    for (p_break in first_b:last_b) {
        tmp_result[p_break] <- ssr_data[beg, p_break] +
            ssr_data[p_break + 1, end]
    }

    tmp_ssr <- min(tmp_result[first_b:last_b])
    tmp_break <- (first_b - 1) + which.min(tmp_result[first_b:last_b])

    return(list(ssr = tmp_ssr, break_point = tmp_break))
}
