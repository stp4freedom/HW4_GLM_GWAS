
if (!file.exists("GAPIT_Tutorial_Data.zip"))
{
  download.file("http://zzlab.net/GAPIT/GAPIT_Tutorial_Data.zip", destfile = "GAPIT_Tutorial_Data.zip")
  unzip("GAPIT_Tutorial_Data.zip")
}
download.file("http://zzlab.net/GAPIT/data/CROP545_Covariates.txt", destfile = "CROPS545_Covariates.txt")
download.file("http://zzlab.net/GAPIT/data/CROP545_Phenotype.txt", destfile = "CROPS545_Phenotype.txt")
# Import the GAPIT demo data genotypes
gt_scan <- data.frame(read.table("GAPIT_Tutorial_Data/mdp_numeric.txt", header = T, stringsAsFactors = F, sep = "\t", nrows = 1))
classes <- sapply(gt_scan, class)
genotypes <- data.frame(read.table("GAPIT_Tutorial_Data/mdp_numeric.txt", header = T, row.names = 1, colClasses = classes, stringsAsFactors = F, sep = "\t"))

GM <- read.table("GAPIT_Tutorial_Data/mdp_SNP_information.txt", header = T, stringsAsFactors = F, sep = "\t")
CV <- read.table("CROPS545_Covariates.txt", header = T, stringsAsFactors = F, sep = "\t")
phenotypes <- read.table("CROPS545_Phenotype.txt", header = T, stringsAsFactors = F, sep = "\t")

GWAStest<- function(phenotypes=NULL, genotypes=NULL, Cov=NULL, GM=NULL, PCA.M=3, QTN.position=NULL, cutoff=NULL){
  print("GWAStest Starting")
  ###check and copy input data
  GD=genotypes
  PCA=prcomp(GD)
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=phenotypes
  P=matrix(NA,1,m)
  for (i in 1:m){
    x=GD[,i]
    # CHECK FOR GENOTYPE DISTRIBUTION
    if(max(x)==min(x)){p=1
    # IF MAX NOT EQUAL TO MIN
    }else{
      # CHECK FOR DEPENDENCE
      # IF NO CV INPUT
      if (is.null(CV)){X=cbind(1, PCA$x[,1:PCA.M],x)
      # IF THERE IS CV INPUT
      }else{fD <- fixDependence(PCA$x[,1:PCA.M], as.matrix(CV))
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
    }# END FOR MAX NOT EQUAL TO MIN
    P[i]=p[length(p)] # COLLECT P-VALUES
  } #end of looping for markers

  P.value=P
  order.SNP=order(P.value)

  #If cutoff is default, uses Bonferroni; else uses -log(value)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]

  ###power.fdr.type1error calculation
  power.fdr.type1error=NULL
  if (!is.null(QTN.position)){
    power.fdr.type1error=power.fdr(order.SNP, QTN.position)
  }
  ###
  zeros=P==0
  P[zeros]=1e-20
  P=data.frame(P)
  ###Generate Manhattan plot
  #color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
  manhattan_plot <- function(marker_map, pvals,cutoff,QTN_index = c(), trait = "unknown")
  {
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
  mp=manhattan_plot(GM,P,cutoff.final)
  knit_print(mp)
  ###Generate QQ plot
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
  qq=qq_plot(GM,P)
  knit_print(qq)
  #myGLM.results = data.frame((c(row.names(t(GD))[-1],GM[,2],GM[,3],P.value,rank(P.value))),nrow=n,ncol=5)
  #colnames(myGLM.results)=c("SNP ","Chromosome ","Position ","P value ","Order ")
  #write.csv(myGLM.results,file="myGLM.results.csv",row.names=FALSE)
  print("GWAStest ran successfully finished!")
  #return(list(P.value=P.value, cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP.P=P.value[sig.SNP], order.SNP=order.SNP, power.fdr.type1error=power.fdr.type1error))
  return(list(cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP.P=P.value[sig.SNP], order.SNP=order.SNP, power.fdr.type1error=power.fdr.type1error))
  }
phenotypes_n_1=phenotypes[,2]
GWAStest(phenotypes_n_1,genotypes,CV,GM,PCA.M = 3)

#type I error
power.fdr <- function(order.SNP, QTN.position) {
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
power.fdr(order.SNP)
