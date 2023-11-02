#***************************************************************************************************
# Aim of this function is extract the p-value matrix from gene-gene causal network
#***************************************************************************************************
# @ - result_list: samples for training the prediciton model
# @ - pcSimple.fit: samples for prediction
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
#***************************************************************************************************
GetPvalueMatrix_ForAdjustment<-function(result_list,pcSimple.fit,g,p)
{
  
  pValue_gg = result_list[[2]]
  
  adj_gg = as(result_list[[1]],"matrix") #update by 9/6/2021.
  pValue_gg[which(adj_gg==0 & pValue_gg < 0.01, arr.ind = T)] = 1 #update by 9/6/2021. Here 0.01 is the alpha from PC-algorithm.
  
  pValue_gp = 2*pnorm(-abs(pcSimple.fit$zMin))
  pValue_rho = matrix(nrow= g+p, ncol=g+p)
  pValue_rho [1:g,1:g] = pValue_gg
  pValue_rho [,g+p] <- c( pValue_gp,-1)
  pValue_rho [g+p,] <- c( pValue_gp,-1)
  
  diag(pValue_rho) <- 0
  pValue_rho[pValue_rho==-1] <- 1
  return(pValue_rho)
}



#***************************************************************************************************
# Aim of this function is adjust the covariance matrix by p-value iterately
#***************************************************************************************************
# @ - iterateNum: the total number of iterations for p-value adjustment
# @ - thresh: threshold for comparing the different between two adjacent iterations
# @ - cov_imputedExpr: covariance matrix based on imputed gene expression data
# @ - pValue_rho: the return value from 'GetPvalueMatrix_ForAdjustment' function
#***************************************************************************************************
Pvalue_adjustment_Iterative<-function(iterateNum = 20,thresh = 0.0001,cov_imputedExpr,pValue_rho)
{
  
  iterate_error2 = array(0, iterateNum)
  cov_imputedExpri2 = cov_imputedExpr
  
  for(i in 1:iterateNum)
  {
    GLASSO = glassoFast(cov_imputedExpri2, rho=pValue_rho)
    cov_imputedExpri2 = GLASSO$w
    iterate_error2[i] = iterate_error2[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
      log(nrow(dat2))/2
    if((i>1)&&(abs(iterate_error2[i]-iterate_error2[i-1])< thresh))
    {
      break
    }
    
  }
  
  cov_imputedExpr5 = cov_imputedExpri2
  return(cov_imputedExpr5)
}