#' Detect linear dependencies of one matrix on another
#'
#' @param X1 A matrix.
#' @param X2 A matrix, the columns of which may be partially linearly dependent on the columns of X1.
#' @param tol The tolerance to use when assessing linear dependence.
#' @param The tolerance to use when assessing linear dependence.
#' @param rank.def If the degree of rank deficiency in X2, given X1, is known, then it can be supplied here, and tol is then ignored. Unused unless positive and not greater than the number of columns in X2.
#' @param strict if TRUE then only columns individually dependent on X1 are detected, if FALSE then enough columns to make the reduced X2 full rank and independent of X1 are detected.
#' @return A vector of the columns of X2 which are linearly dependent on columns of X1 (or which need to be deleted to acheive independence and full rank if strict==FALSE). NULL if the two matrices are independent.
fix_Dep<-function (X1, X2, tol = .Machine$double.eps^0.5, rank.def = 0,
          strict = FALSE)
{
  qr1 <- qr(X1, LAPACK = TRUE)
  R11 <- abs(qr.R(qr1)[1, 1])
  r <- ncol(X1)
  n <- nrow(X1)
  if (strict) {
    QtX2 <- qr.qty(qr1, X2)
    QtX2[-(1:r), ] <- 0
    mdiff <- colMeans(abs(X2 - qr.qy(qr1, QtX2)))
    if (rank.def > 0)
      ind <- (1:ncol(X2))[rank(mdiff) <= rank.def]
    else ind <- (1:ncol(X2))[mdiff < R11 * tol]
    if (length(ind) < 1)
      ind <- NULL
  }
  else {
    QtX2 <- qr.qty(qr1, X2)[(r + 1):n, ]
    qr2 <- qr(QtX2, LAPACK = TRUE)
    R <- qr.R(qr2)
    r0 <- r <- nrow(R)
    if (rank.def > 0 && rank.def <= nrow(R))
      r0 <- r - rank.def
    else while (r0 > 0 && mean(abs(R[r0:r, r0:r])) < R11 *
                tol) r0 <- r0 - 1
    r0 <- r0 + 1
    if (r0 > r)
      return(NULL)
    else ind <- qr2$pivot[r0:r]
  }
  ind
}
