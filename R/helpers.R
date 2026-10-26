asterisk <- function(x, k = c(0.10, 0.05, 0.01, 0.001)) {
  stopifnot(length(k) > 0, all(diff(k) < 0), !anyNA(k))

  n <- length(k)
  i <- 0:n

  symbols <- if (n <= 3) {
    strrep("*", i)
  } else if (n == 4) {
    c("", ".", strrep("*", seq_len(n - 1)))
  } else {
    paste0(strrep("*", i %/% 2), ifelse(i %% 2 == 1, ".", ""))
  }

  level <- n - findInterval(x, rev(k), left.open = TRUE)
  symbols[level + 1]
}

testcols <- function(x) {
  m <- crossprod(as.matrix(x))
  ee <- eigen(m)
  evecs <- split(zapsmall(ee$vectors), col(ee$vectors))
  mapply(function(val, vec) {
    if (val != 0) NULL else which(vec != 0)
  }, zapsmall(ee$values), evecs)
}
