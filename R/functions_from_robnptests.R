#' @export
win_mean <- function (x, gamma = 0.2, na.rm = FALSE) {
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, 
                            null.ok = FALSE)
  checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, 
                           upper = 0.5, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  }
  else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }
  if (identical(gamma, 0)) {
    return(mean(x))
  }
  n <- length(x)
  r <- floor(gamma * n)
  x.sort <- sort(x)
  x.sort[1:r] <- x.sort[r + 1]
  x.sort[(n - r + 1):n] <- x.sort[n - r]
  return(mean(x.sort))
}

#' @export
win_var <- function (x, gamma = 0, na.rm = FALSE) {
  stopifnot(`'x' is missing.` = !missing(x))
  checkmate::assert_numeric(x, min.len = 2, finite = TRUE, 
                            all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_number(gamma, lower = 0, upper = 0.5, na.ok = FALSE, 
                           finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  if (!na.rm & any(is.na(x))) {
    return(list(var = NA_real_, h = NA_real_))
  }
  else if (na.rm) {
    x <- as.vector(stats::na.omit(x))
  }
  if (length(unique(x)) == 1) {
    stop("All values in '", deparse(substitute(x)), 
         "' are equal. The scale estimate is '0' and the test statistic cannot be computed.")
  }
  n <- length(x)
  r <- floor(gamma * n)
  x.sort <- sort(x)
  x.lower <- x.sort[r + 1]
  x.upper <- x.sort[n - r]
  x.sort[which(x.sort < x.lower)] <- x.lower
  x.sort[which(x.sort > x.upper)] <- x.upper
  res <- 1/(n - 1) * sum((x.sort - mean(x.sort))^2)
  h <- n - 2 * r
  #return(list(var = res, h = h))
  return(res)
}

#' @export
win_sd <- function(x, ...) {
 
  return(sqrt( win_var(x, ...) ))
  
}

#' @export
trim_mean <- function (x, gamma = 0.2, na.rm = FALSE) {
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, 
                            null.ok = FALSE)
  checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, 
                           upper = 0.5, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  return(mean(x, trim = gamma, na.rm = na.rm))
}