#' Function to create a QQ plot
#'
#' @param GM genetic map of data with chr and position of each SNP
#' @param pvals pvals from gwas results for each SNP
#' @param QTN_index posistion of QTN if applicable
#' @param trait character value for trait name
#' @return QQ plot
qq_plot <- function(GM, pvals, QTN_index= c(), trait = "unknown")
{
  GM$pvals <- as.vector(-log10(t(pvals))) # Add pvalues to the data.frame and log10 transform

  exp_pval_dist <- -log10(runif(nrow(GM), 0, 1)) # Sample random p-values from a uniform distribution between 0 and 1

  qq_plot <- ggplot(GM, aes(x = sort(exp_pval_dist),
                                    y = sort(GM$pvals))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = paste("Q-Q Plot for trait:", as.character(trait)),
         x = "-log10(p-values) Expected",
         y = "-log10(p-values) Observed")

  return(qq_plot)
}
