#' Outlier detection function
#' @param x numeric vector
#' @param z standardized threshold
#' @return vector index of outliers
#' @export
#' @examples
#' out(runif(10), 0.9)
out <- function (x, z = 0.7)
    which(abs(scale(x, scale = TRUE, center = TRUE)) >= z)

