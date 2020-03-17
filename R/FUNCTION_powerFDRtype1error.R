#' Function to caculate power, FDR and type-1 error
#'
#' @param order.SNP list of SNPs order by ascending  p-value
#' @param QTN.position position of QTN if known
#' @return list of power, FDR and type-1 error
power.fdr <- function(order.SNP, QTN.position=NULL) {
  pwr <- c()
  fdr <- c()
  t1error <- c()
  NQTN <- length(QTN.position)
  nsnp <- length(order.SNP)
  for (m in (1: nsnp)) {
    detected <- intersect(order.SNP[1:m], QTN.position)
    falsePositive <- setdiff(order.SNP[1:m], QTN.position)
    pwr <- c(pwr, length(detected)/NQTN)
    fdr <- c(fdr, length(falsePositive)/m)
    t1error <- c(t1error, length(falsePositive)/nsnp)
  }
  return(list(power=pwr, fdr=fdr, type1error=t1error))
}

