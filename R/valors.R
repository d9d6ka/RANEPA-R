#' @export
valors <- function(model, lambda, k) {
    m_vc <- coi_end_vc[[model]][[k]]
    if (0 < lambda & lambda <= 0.15)
        v_vc <- m_vc[, 1, drop = FALSE]
    else if (0.15 < lambda & lambda <= 0.25)
        v_vc <- m_vc[, 2, drop = FALSE]
    else if (0.25 < lambda & lambda <= 0.35)
        v_vc <- m_vc[, 3, drop = FALSE]
    else if (0.35 < lambda & lambda <= 0.45)
        v_vc <- m_vc[, 4, drop = FALSE]
    else if (0.45 < lambda & lambda <= 0.55)
        v_vc <- m_vc[, 5, drop = FALSE]
    else if (0.55 < lambda & lambda <= 0.65)
        v_vc <- m_vc[, 6, drop = FALSE]
    else if (0.65 < lambda & lambda <= 0.75)
        v_vc <- m_vc[, 7, drop = FALSE]
    else if (0.75 < lambda & lambda <= 0.85)
        v_vc <- m_vc[, 8, drop = FALSE]
    else if (0.85 < lambda & lambda < 1)
        v_vc <- m_vc[, 9, drop = FALSE]
    else
        stop("ERROR! Try to specify the value of lambda again")

    return(v_vc)
}
