#' Outlier detection function
#' @param x numeric vector
#' @param z standardized threshold
#' @return vector index of outliers
out <- function (x, z = 0.7)
    which(abs(scale(x, scale = TRUE, center = TRUE)) >= z)
