#' Function to create a manhattan plot
#'
#' @param GM genetic map of data with chr and position of each SNP
#' @param pvals pvals from gwas results for each SNP
#' @param cutoff  If cutoff is default, uses Bonferroni; else uses -log(value) of 0.05/number of SNPs
#' @param QTN_index posistion of QTN if applicable
#' @param FP Positions of SNPS that are False positive
#' @param TP Positions of SNPS that are True positive
#' @param trait character value for trait name
#' @return Manhatten plot
manhattan_plot <- function(GM, pvals,cutoff=NULL,QTN_index = c(),FP=NULL,TP=NULL, trait = "unknown",messages=FALSE){
  m=nrow(GM)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  if(messages==TRUE){print(paste0("The final cuttoff for a significant p-value is ",as.character(cutoff.final)))}
  GM$pvals <- -log10(t(pvals)) # Add pvalues to the data.frame and log10 transform
  SG_pval <- -log10(cutoff.final)
  GM$comb_pos <- GM$Chromosome * 1e9 + GM$Position
  if(is.null(FP)){
  manhattan_plot <- ggplot(GM, aes(x = 1:nrow(GM), y = pvals, color = factor(Chromosome))) +
    geom_point() +
    geom_vline(xintercept = QTN_index, color = "red") +
    geom_hline(yintercept= SG_pval,color="green")+
    labs(title = paste("GWAS manhattan plot for trait:", as.character(trait)),
         y = "-log10(p-value)",
         x = "Marker Position",
         color = "Chromosome")

  return(manhattan_plot)
  }else{
    if((!is.null(TP))){
      manhattan_plot <- ggplot(GM, aes(x = 1:nrow(GM), y = pvals, color = factor(Chromosome))) +
        geom_point() +
        geom_vline(xintercept = QTN_index, color = "red") +
        geom_vline(xintercept = FP, color = "blue") +
        geom_vline(xintercept = TP, color = "black") +
        geom_hline(yintercept= SG_pval,color="green")+
        labs(title = paste("GWAS manhattan plot for trait:", as.character(trait)),
             y = "-log10(p-value)",
             x = "Marker Position",
             color = "Chromosome")
      return(manhattan_plot)
    }else{  manhattan_plot <- ggplot(GM, aes(x = 1:nrow(GM), y = pvals, color = factor(Chromosome))) +
      geom_point() +
      geom_vline(xintercept = QTN_index, color = "red") +
      geom_vline(xintercept = FP, color = "blue") +
      geom_hline(yintercept= SG_pval,color="green")+
      labs(title = paste("GWAS manhattan plot for trait:", as.character(trait)),
           y = "-log10(p-value)",
           x = "Marker Position",
           color = "Chromosome")
    return(manhattan_plot)}

  }
}


