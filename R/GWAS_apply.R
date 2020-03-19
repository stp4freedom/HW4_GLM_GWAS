#' GWAS with PCA
#'
#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param Cov numeric data.frame with covariates values
#' @param GM genetic map of data with chr and position of each SNP
#' @param PCA.M number of principal components to use default is 3
#' @param QTN.position posistion of QTN if applicable
#' @param cutoff  If cutoff is default, uses Bonferroni; else uses -log(value) of 0.05/number of SNPs
#' @param plots  if TRUE, function plots PCA graphs, Manhatten Plot and QQ plot
#' @param messages if TRUE, returns messages for the GWAS function
#' @param print if TRUE, results are saved in a CSV
#' @param trait character value for trait name
#' @return Manhatten plot, QQ plot plus p-values, type-1 error and power for every SNP and results in a CSV
GWASapply<- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL, PCA.M=3, QTN.position=NULL, cutoff=NULL,plots=FALSE,messages=FALSE,print=FALSE,trait="unknown"){
  #If messages is true display gwas starting
  if(messages==TRUE){print("GWASapply Starting")}
  apply_start <- proc.time()
  ###check and copy input data
  GD=geno
  PCA=prcomp(GD)
  if(messages==TRUE){print("Principal Components have been calculated Successfully")}
  PCA_results <- summary(PCA)
  if(plots==TRUE){print(kable(round(PCA_results$importance[,1:10], 2)))

  #Variance Explained
  var_exp_plot_data <- data.frame(t(PCA_results$importance)) # I transposed the data to be in long form and coerce to a data.frame for ggplot
  names(var_exp_plot_data) <- c("sdev", "var_exp", "cum_var_exp") # rename the columns (because the original rownames has spaces)
  var_exp_plot_data$pca_index <- 1:nrow(var_exp_plot_data) # Add a new column that is the integer index of the component
  var_exp_plot <- ggplot(data = var_exp_plot_data, aes(x = pca_index, y = var_exp)) +
    geom_line() +
    geom_point() +
    labs(x = "PCA component index", y = "Variance explained", title = "Variance Explained")

  #Cumulative Variance Explained
  cum_var_exp_plot <- ggplot(data = var_exp_plot_data, aes(x = pca_index, y = cum_var_exp)) +
    geom_line() +
    geom_point() +
    labs(x = "PCA component index", y = "Cumulative Variance explained", title = "Cumulative Variance Explained")
  grid.arrange(var_exp_plot, cum_var_exp_plot, nrow = 1, ncol = 2, top = "Variance of Principal Components")
  #normal_print(pcav)
  # Plotting the first three components
  PCA_plot_data <- data.frame(PCA$x)
  pca_comp_plot_12 <-
    ggplot(data = PCA_plot_data, aes(x = PC1, y = PC2)) +
    geom_point()
  pca_comp_plot_13 <-
    ggplot(data = PCA_plot_data, aes(x = PC1, y = PC3)) +
    geom_point()
  pca_comp_plot_23 <-
    ggplot(data = PCA_plot_data, aes(x = PC2, y = PC3)) +
    geom_point()
  grid.arrange(pca_comp_plot_12, pca_comp_plot_13, pca_comp_plot_23, nrow = 2, ncol = 2, top = "Principal Components")
  if(messages==TRUE){print("Principal Components plots have been printed Successfully")}}
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=pheno
  P=apply(GD,2, function(x)

    #x=GD[,1]
    # CHECK FOR GENOTYPE DISTRIBUTION
    if(max(x)==min(x)){p=1
    P=p[length(p)]
    # IF MAX NOT EQUAL TO MIN
    }else{
      # CHECK FOR DEPENDENCE
      # IF NO CV INPUT
      if (is.null(CV)){X=cbind(1, PCA$x[,1:PCA.M],x)
      # IF THERE IS CV INPUT
      }else{fD <- fix_Dep(PCA$x[,1:PCA.M], as.matrix(CV))
      # WITH CV INPUT, CONDITION 1: NO DEPENDENCE
      if (is.null(fD)){X=cbind(1, PCA$x[,1:PCA.M],CV, x)

      # WITH CV INPUT, CONDITION 2: WITH DEPENDENCE
      }else {X=cbind(1, PCA$x[,1:PCA.M],x)}
      }# END FOR DEPENDECE
      # SOLVE THE LINEAR REGRESSION MATRIX:
      #  X=as.matrix(X)
      LHS=t(X)%*%X
      C=solve(LHS)
      RHS=t(X)%*%y
      b=C%*%RHS
      yb=X%*%b
      e=y-yb
      n=length(y)
      ve=sum(e^2)/(n-1)
      vt=C*ve
      t=b/sqrt(diag(vt))
      p=2*(1-pt(abs(t),n-2))
      P=p[length(p)]
    }# END FOR MAX NOT EQUAL TO MIN
    # COLLECT P-VALUES
  ) #end of appply function for markers
  P=t(matrix(P))
  P.value=P
  order.SNP=order(P.value)
  #If cutoff is default, uses Bonferroni; else uses -log(value)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  if(messages==TRUE){print(paste0("The final cuttoff for a significant p-value is ",as.character(cutoff.final)))}
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]
  lsnp=length(sig.SNP)
  if(messages==TRUE){print(paste0(as.character(lsnp), " Significant SNPs were found"))}

  ###
  zeros=P==0
  P[zeros]=1e-20
  P=data.frame(P)
  ###Generate Manhattan plot
  if (is.null(QTN.position)){

    if(plots==TRUE){mp=manhattan_plot(GM,P,cutoff.final,trait=as.character(trait))
  knit_print(mp)
  if(messages==TRUE){print("Manhattan Plot Printed without QTN")}}

  # WITH CV INPUT, CONDITION 2: WITH DEPENDENCE
  }else {
    ###power.fdr.type1error calculation
    power.fdr.type1error=NULL
    if (!is.null(QTN.position)){
      power.fdr.type1error=power.fdr(P.value, QTN.position,cutoff.final)
    }
    detected=intersect(sig.SNP,QTN.position)
    if(length(detected)==0){detected=NULL}
    falsePositive=setdiff(sig.SNP, QTN.position)
    if(messages==TRUE){print(paste0("QTN's detected with ",as.character(length(detected))," True Positives and ",as.character(length(falsePositive))," False Positives."))}
    if(plots==TRUE){ mp=manhattan_plot(GM,P,cutoff.final,QTN.position,falsePositive,detected,trait=as.character(trait))
    knit_print(mp)
    if(messages==TRUE){print("Manhattan Plot Printed with QTNs and False and True Positives")}}

  }

  ###Generate QQ plot
  if(plots==TRUE){qq=qq_plot(GM,P,trait=as.character(trait))
  knit_print(qq)
  if(messages==TRUE){print("QQ-Plot Printed")}}
  P.value_df=data.frame(P.value)
  if(print==TRUE){GLM.results = data.frame(GM,P.value_df,rank(P.value_df))
  colnames(GLM.results)=c("SNP ","Chromosome ","Position ","P value","Order")
  fwrite(GLM.results,file="GLM.results.csv",row.names=FALSE)}
  if(messages==TRUE){print("GWASapply results have printed")}
  #Returns Values that include false and true positives if QTN are provided
  if (!is.null(QTN.position)){
    GWAS.Results=list(True.Positive=detected,False.Positive=falsePositive,P.value.res=P.value, cutoff.final.res=cutoff.final, sig.SNP.res=sig.SNP, sig.SNP.P.res=P.value[sig.SNP], order.SNP.res=order.SNP, power.fdr.type1error.res=power.fdr.type1error)
  }else{
    GWAS.Results=list(P.value.res=P.value, cutoff.final.res=cutoff.final, sig.SNP.res=sig.SNP, sig.SNP.P.res=P.value[sig.SNP], order.SNP.res=order.SNP)
  }
  apply_end <- proc.time()
  apply_elapsed <- apply_end[3] - apply_start[3]
  if(messages==TRUE){print(paste0("GWASapply ran successfully and took ",as.character(round(apply_elapsed[[1]],2))," seconds"))}
  return(GWAS.Results)
}

