#' Function to caculate power, FDR and type-1 error, False Positives, and True Positives
#'
#' @param P list of SNPs order by ascending  p-value
#' @param QTN.position position of QTN if known
#' @return list of power, FDR, type-1 error, False Positives, and True Positives
power.fdr <- function(P.value, QTN.position=NULL,cutoff=NULL) {
  pwr <- c()
  fdr <- c()
  t1error <- c()
  NQTN <- length(QTN.position)
  nsnp <- length(P.value)
  order.SNP=order(P.value)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]
  TP=intersect(sig.SNP,QTN.position)
  FP=setdiff(sig.SNP, QTN.position)
  for (m in (1: nsnp)) {
    detected <- intersect(order.SNP[1:m], QTN.position)
    falsePositive <- setdiff(order.SNP[1:m], QTN.position)
    pwr <- c(pwr, length(detected)/NQTN)
    fdr <- c(fdr, length(falsePositive)/m)
    t1error <- c(t1error, length(falsePositive)/nsnp)
  }
  return(list(power=pwr, fdr=fdr, type1error=t1error,FP.fdr.power=FP,TP.fdr.power=TP))
}

