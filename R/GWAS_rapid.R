GWASapply_rapid<- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL, PCA.M=3,cutoff=NULL){
  GD=geno
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=pheno
  PCA=prcomp(GD)
  P=apply(GD,2, function(x)
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
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]
  lsnp=length(sig.SNP)
  ###
  zeros=P==0
  P[zeros]=1e-20
  P=data.frame(P)
  GWAS.Results=list(P.value.res=P.value, cutoff.final.res=cutoff.final, sig.SNP.res=sig.SNP, sig.SNP.P.res=P.value[sig.SNP], order.SNP.res=order.SNP)
  return(GWAS.Results)
}
