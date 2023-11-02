#***************************************************************************************************
# Aim of this function is to impute the gene expression data based on the information of genotypes
#***************************************************************************************************
# @ - samples_data: samples for prediction
# @ - cross_validation_data: samples for training the model
#***************************************************************************************************
Imupte_GeneExpressionData<-function(samples_data,cross_validation_data)
{
  # allocate snps data and gene data randomly
  snps_data = cross_validation_data[,snp_index]  ##the 800 subjects (GTEx)
  gene_data = cross_validation_data[,gene_index] ##the 800 subjects (GTEx)
  snps_data2 = samples_data[,snp_index] #use for prediction
  
  
  samples_data_iptg  = samples_for_impute
  
 
  for (k in gene_index)
  {
    
    gene_k = cross_validation_data[,k] # obtain the gene data in turn 
    
    cvfit =cv.glmnet(snps_data,gene_k,nfolds = cv.glmnet_nfolds)
    
    coef(cvfit, s = "lambda.min")
    
    
    #choosing the best lambda according to cross-validation by using  cv.glmnet
    pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = "lambda.min")
    if (var(pred_genes)==0) {
      
      cvm_min_nonzero = min( cvfit$cvm[cvfit$nzero>0] )  #extract the best cv error at which there is at least one non-zero coeff. 
      lamb = cvfit$lambda[cvfit$cvm==cvm_min_nonzero]
      pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = lamb)
    }
    
  
    samples_data_iptg[,k] = pred_genes[,1:1]  #use the 1st column data temporarily
    
  
  }
  
  return(samples_data_iptg)
  
}




#***************************************************************************************************
# Aim of this function is to impute the data of Trait1 based on genotype information
#***************************************************************************************************
# @ - training_data: samples for training the prediciton model
# @ - pred_data: samples for prediction
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index: index of snps in the simulated causal network
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
# @ - cv.glmnet_nfolds: number of folds for cross validation
#***************************************************************************************************
Imupte_T1Data<-function(training_data,pred_data,gene_index,snp_index,s,g,p,cv.glmnet_nfolds=5)
{
  snps_data = training_data[,snp_index]
  T1_index = s+g+1
  T1_data = training_data[,T1_index]
  cvfit =cv.glmnet(snps_data,T1_data,nfolds = cv.glmnet_nfolds)

  cf<-coef(cvfit, s = "lambda.min")
  i<-which(cvfit$lambda == cvfit$lambda.min)
  e<-cvfit$cvm[i]
  rsq <-1-e/var(T1_data)
  
  snps_data_pred = pred_data[,snp_index]
  
  #choosing the best lambda according to cross-validation by using cv.glmnet
  pred_T1 = predict(cvfit, newx = snps_data_pred, s = "lambda.min")
  
  T2_index = s+g+2
  dt1 = pred_data[,gene_index]
  colnames(dt1) = seq(1,ncol(dt1))
  
  dt2 = pred_T1
  dt2 = as.data.table(dt2)
  colnames(dt2) = "T1"
  
  
  dt3 = pred_data[,T2_index]
  dt3 = as.data.table(dt3)
  colnames(dt3) = "T2"

  dat = cbind(dt1,dt2,dt3)
  res_list = list()
  res_list[["dat"]] = dat
  res_list[["rsq"]] = rsq
  return(res_list)
  
}




#***************************************************************************************************
# Aim of this function is to get the correlation of imputed gene expression data and raw gene expression data
#***************************************************************************************************
# @ - gene_index: index of genes in the simulated causal network
# @ - samples_data: samples which include raw gene expression data
# @ - samples_data_iptg: samples which include imputed gene expression data
# @ - sample_num: the number of samples in total
#***************************************************************************************************
Caculate_Cor_ImputeGeneAndRawGene<-function(gene_index,samples_data,samples_data_iptg,sample_num)
{
  dt_genedata = data.frame()
  row_count = 1
  
  for (g_count in gene_index) # iterate over each gene according to gene_index.
  {

    v1 = samples_data[,g_count]
    v2 = samples_data_iptg[,g_count]
    
    # results in comparsion
    MSE_value = MSE(v1,v2)/sample_num 
    cor_value = cor(v1,v2)  # in some cases, cor value is NA, because the standard deviation is zero
    cosine_value = cosine(v1,v2) 
    
    dt_genedata[row_count,1] = paste0("G",g_count)
    dt_genedata[row_count,2] = MSE_value
    dt_genedata[row_count,3] = cor_value
    dt_genedata[row_count,4] = cosine_value
    colnames(dt_genedata)[1] = "Gene_num"
    colnames(dt_genedata)[2] = "MSE_value"
    colnames(dt_genedata)[3] = "cor_value"
    colnames(dt_genedata)[4] = "cosine_value"
    row_count = row_count+1
  }
  
  return(dt_genedata) #print gene imputed data and simulated data
  
}




#***************************************************************************************************
# Aim of this function is learn the gene-pheno network with no-hidden variables
#***************************************************************************************************
# @ - gene_index: index of genes in the simulated causal network
# @ - samples_data_iptg: samples which include imputed gene expression data
#***************************************************************************************************
GeneToPhenotype_Network_nohidden<-function(samples_data_iptg,gene_index)
{
  
  var_thres = 0.01
  column3 = s+g+1
  column4 = s+g+p
  
  pheno_data = as.matrix(samples_data_iptg[1:sample_num_predict_imputedgene, column3:column4])
  gene_data = samples_data_iptg[1:sample_num_predict_imputedgene,gene_index] # extract the gene data according to gene_index
  
  dfnew_rankNorm = cbind(pheno_data,gene_data) 
  colnames(dfnew_rankNorm)[1] <-"outcome"
  
  
  library(coop) 
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat) 
  
  
  #********************************************************
  # Apply PC-simple algorithm around the response variable on NFBC
  #**************************************************************
  library(pcalg)
  
  ## Load  data
  n <- nrow (dfnew_rankNorm)
  V <- colnames(dfnew_rankNorm) # labels aka node names
  
  t1 = proc.time()
  
  pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1], 
                                    dm=dfnew_rankNorm[,-1] , 
                                    corMat = precompute_corMat, 
                                    alpha=0.001, 
                                    max_ord=3, 
                                    corMethod = "standard",
                                    verbose = TRUE, directed = TRUE)
  proc.time()-t1
  return(pcSimple.fit)
}


#***************************************************************************************************
# Aim of this function is learn the gene-pheno network with multiple traits in the simulated causal network
#***************************************************************************************************
# @ - samples_data_iptg: samples which include imputed gene expression data
#***************************************************************************************************
GeneToPhenotype_Network_MultiTraits<-function(samples_data_iptg)
{
  
  var_thres = 0.01
  column3 = s+g+1
  p=1 #Added by fymn. 2023.9.25. Force the number of pheno equal to 1, for each time we learn gene-pheno network
  column4 = s+g+p
  
  pheno_data = as.matrix(samples_data_iptg[1:sample_num_predict_imputedgene, column3:column4])
  
  
  gene_data = samples_data_iptg[1:sample_num_predict_imputedgene,gene_index] # extract the gene data according to gene_index
  
  
  dfnew_rankNorm = cbind(pheno_data,gene_data) #鐠嬪啯宕瞘ene閸滃henotype閻ㄥ嫭鏆熼幑顕嗙礉閹跺henotype閺佺増宓侀弨鎯ф躬閸撳秹�??
  colnames(dfnew_rankNorm)[1] <-"outcome"
  
  
  library(coop) 
  
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  
  
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat) 
  
  
  #********************************************************
  # Apply PC-simple algorithm around the response variable on NFBC
  #**************************************************************
  library(pcalg)

  n <- nrow (dfnew_rankNorm)
  V <- colnames(dfnew_rankNorm) # labels aka node names
  
  t1 = proc.time()
    
  pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1], 
                                    dm=dfnew_rankNorm[,-1] , 
                                    corMat = precompute_corMat, 
                                    alpha=0.001, 
                                    max_ord=3, 
                                    corMethod = "standard",
                                    verbose = TRUE, directed = TRUE)
  proc.time()-t1
  
  return(pcSimple.fit)
}




#***************************************************************************************************
# Aim of this function is learn the gene-pheno network by using imputed trait1 data and raw trait2 data
#***************************************************************************************************
# @ - ImputedT1: samples for learning gene-pheno network with imputed trait1 data
#***************************************************************************************************
GeneToPhenotype_Network_geneT2<-function(ImputedT1)
{
  index_T1 = ncol(ImputedT1)-1
  index_T2 = ncol(ImputedT1)
  
  pheno_data = ImputedT1[,index_T2:index_T2]
  gene_data = ImputedT1[,1:index_T1]
  
  dfnew_rankNorm = cbind(pheno_data,gene_data) 
  colnames(dfnew_rankNorm)[1] <-"outcome"
  colnames(dfnew_rankNorm)[2:ncol(dfnew_rankNorm)] = seq(1,(ncol(dfnew_rankNorm)-1))
  
  dfnew_rankNorm = as.matrix(dfnew_rankNorm)
  
  library(coop) 
  
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat) 
  
  
  #********************************************************
  # Apply PC-simple algorithm around the response variable on NFBC
  #**************************************************************
  library(pcalg)

  n <- nrow (dfnew_rankNorm)
  V <- colnames(dfnew_rankNorm) # labels aka node names
  
  t1 = proc.time()
  ## estimate local causal network structure around the response variable 

  pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1], 
                                    dm=dfnew_rankNorm[,-1] , 
                                    corMat = precompute_corMat, 
                                    alpha=0.001, 
                                    max_ord=3, 
                                    corMethod = "standard",
                                    verbose = TRUE, directed = TRUE)
  
  
  proc.time()-t1
  
  return(pcSimple.fit)
}
