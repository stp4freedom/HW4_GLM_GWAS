qq_plot <- function(marker_map, pvals, QTN_index= c(), trait = "unknown")
{
  marker_map$pvals <- -log10(t(pvals)) # Add pvalues to the data.frame and log10 transform

  exp_pval_dist <- -log10(runif(nrow(marker_map), 0, 1)) # Sample random p-values from a uniform distribution between 0 and 1

  qq_plot <- ggplot(marker_map, aes(x = sort(exp_pval_dist),
                                    y = sort(marker_map$pvals))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "Q-Q Plot for GWAStest",
         x = "-log10(p-values) Expected",
         y = "-log10(p-values) Observed")

  return(qq_plot)
}
