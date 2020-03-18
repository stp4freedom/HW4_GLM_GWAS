manhattan_plot <- function(marker_map, pvals,cutoff,QTN_index = c(), trait = "unknown"){
  marker_map$pvals <- -log10(t(pvals)) # Add pvalues to the data.frame and log10 transform
  SG_pval <- -log10(cutoff)
  marker_map$comb_pos <- marker_map$Chromosome * 1e9 + marker_map$Position
  manhattan_plot <- ggplot(marker_map, aes(x = 1:nrow(marker_map), y = pvals, color = factor(Chromosome))) +
    geom_point() +
    geom_vline(xintercept = QTN_index, color = "red") +
    geom_hline(yintercept= SG_pval,color="green")+
    labs(title = paste("GWAS manhattan plot for GWAStest:", trait),
         y = "-log10(p-value)",
         x = "Marker Position",
         color = "Chromosome")

  return(manhattan_plot)
}
