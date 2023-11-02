
### Please ensure that you have modified path to your local path in order to load the according functions
source("/exeh_4/yaning_feng/04_Simulation/R_Package/Simulation_Extension_Functions/Simulation_ExtensionFunction_Mediator_functions.R")



#***************************************************************************************************
# Aim of this function is to organize the results of IDA and jointIDA
#***************************************************************************************************
# @ - dt_IDA: the results of estimated IDA and jointIDA based on our causal graph
#***************************************************************************************************
GetMoreDetail_OverallPerformance<-function(dt_IDA)
{
  #obtain minimum jointIDA, result in est_jointIDA
  est_jointIDA = NULL
  ind1 = which(colnames(dt_IDA)=="true_jointIDA") # obtain the column No of true_jointIDA
  ind2 = which(colnames(dt_IDA)=="Rank_ABS_true_IDA") # obtain the column No of Rank_ABS_true_IDA
  jointIDA_cols = as.data.frame( dt_IDA[,(ind1+1):(ind2-1)], nrow= nrow(dt_IDA)) # all the alternative value of tureIDA
  jointIDA_cols[jointIDA_cols=="NA"] <- 99999
  for (j in 1:nrow(jointIDA_cols) ) {
    rowi = as.numeric(as.character(jointIDA_cols[j,]))
    min_index = which.min(abs(rowi))
    est_jointIDA[j] = jointIDA_cols[j,min_index]
  }
  est_jointIDA = as.numeric(as.character(est_jointIDA))
  dt_IDA = data.frame(dt_IDA, est_jointIDA)
  
  
  #obtain minimum IDA, result in Estimate_IDA_final
  Estimate_IDA_final = NULL
  ind1 = which(colnames(dt_IDA)=="Estimate_IDA")
  ind2 = which(colnames(dt_IDA)=="true_jointIDA")
  IDA_cols = as.data.frame( dt_IDA[,ind1:(ind2-1)], nrow= nrow(dt_IDA)) 
  IDA_cols[IDA_cols=="NA"] <- 99999
  for (j in 1:nrow(IDA_cols) ) {
    rowi = as.numeric(as.character(IDA_cols[j,]))
    min_index = which.min(abs(rowi))
    Estimate_IDA_final[j] = IDA_cols[j,min_index]
  }
  Estimate_IDA_final = as.numeric(as.character(Estimate_IDA_final))	
  dt_IDA = data.frame(dt_IDA, Estimate_IDA_final)
  return(dt_IDA)
}




#***************************************************************************************************
# Aim of this function is to caculate the correlation and RMSE based on the value of IDA and jointIDA
#***************************************************************************************************
# @ - dt_IDA: the results of estimated IDA and jointIDA based on our causal graph
# @ - adj: adjacency matrix of true graph
# @ - g1: the simulated true causal graph
# @ - g2: the estimated causal graph based on our framework
# @ - g3: a random causal graph which is used for comparison
# @ - gene_index: index of genes in the simulated causal network
#***************************************************************************************************
Get_related_variable<-function(dt_IDA,adj,g1,g2,g3,gene_index)
{
  #save correlation of IDA
  Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final,  use='complete.obs' ) 
  Cor_TureIDA_EstIDA_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final, use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson	= cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs')
  Cor_TureIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs', method="spearman")
  Cor_TureIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof, method="spearman")
  
  
  #save correlation of nonzero IDA
  ind_IDA =  which(abs(dt_IDA$true_IDA)>1e-10)
  Cor_TureIDA_EstIDA_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA_final[ind_IDA],  use='complete.obs' )
  Cor_TureIDA_EstIDA_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA_final[ind_IDA], use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs' )
  Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs', method="spearman" )
  Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA])
  Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA], method="spearman")
  
  #save RMSE
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final)
  RMSE_TrueIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueIDA_rowGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  
  #____________________________
  # correlation of joint IDA
  #____________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  # save correlation of joint IDA
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  Cor_TureJointIDA_EstJointIDA_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA, method="spearman")
  
  # modified by fyn 14/2/2020
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof,  method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureJointIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof,  method="spearman")
  
  
  #________________________________________________
  # correlation of joint IDA which are non-zero
  #________________________________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  ind_nonzero =  which(abs(as.numeric(dt_IDA$true_jointIDA))>1e-10)
  
  Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero])
  Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero], method="spearman")
  
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero])
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero],method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero])	
  Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero],  method="spearman")
  
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  RMSE_TrueJointIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueJointIDA_rowGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  # proportion of gene-gene connection
  gene_mat = adj[gene_index,gene_index]
  sum(gene_mat!=0)/nrow(gene_mat)^2
  
  proportion_GeneToGene_Net = sum(gene_mat!=0)/nrow(gene_mat)^2
  
  # correlation of regression coef derived from real expr vs imputed expr
  Cor_rowGene_ImputeGene_Rcof = cor(dt_IDA$Gene_rowdata_regression_Cof, dt_IDA$Gene_imputedata_regression_Cof)
  
  # shd_ture_estGrapph = shd(g1,g2)
  # shd_ture_randomGraph = shd(g1,g3)
  shd_ture_estGrapph = 0
  shd_ture_randomGraph = 0
  
  # correlation and MSE
  related_variable = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Spearman,Cor_TureIDA_ImputeGene_Rcof_Pearson,Cor_TureIDA_ImputeGene_Rcof_Spearman,	Cor_TureIDA_rowGene_Rcof_Pearson,Cor_TureIDA_rowGene_Rcof_Spearman,Cor_TureIDA_EstIDA_Pearson_Nonzero,	Cor_TureIDA_EstIDA_Spearman_Nonzero,	Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero,Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueIDA_EstIDA,	RMSE_TrueIDA_ImputeGene_Rcof,	RMSE_TrueIDA_rowGene_Rcof,	Cor_TureJointIDA_EstJointIDA_Pearson,	Cor_TureJointIDA_EstJointIDA_Spearman,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman,	Cor_TureJointIDA_rowGene_Rcof_Pearson,	Cor_TureJointIDA_rowGene_Rcof_Spearman,	Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero,	Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueJointIDA_EstJointIDA ,	RMSE_TrueJointIDA_ImputeGene_Rcof ,	RMSE_TrueJointIDA_rowGene_Rcof ,	proportion_GeneToGene_Net,	Cor_rowGene_ImputeGene_Rcof,	shd_ture_estGrapph,	shd_ture_randomGraph)
  related_variable_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Spearman","Cor_TureIDA_ImputeGene_Rcof_Pearson","Cor_TureIDA_ImputeGene_Rcof_Spearman","Cor_TureIDA_rowGene_Rcof_Pearson","Cor_TureIDA_rowGene_Rcof_Spearman","Cor_TureIDA_EstIDA_Pearson_Nonzero","Cor_TureIDA_EstIDA_Spearman_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_ImputeGene_Rcof","RMSE_TrueIDA_rowGene_Rcof","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Spearman","Cor_TureJointIDA_ImputeGene_Rcof_Pearson","Cor_TureJointIDA_ImputeGene_Rcof_Spearman","Cor_TureJointIDA_rowGene_Rcof_Pearson","Cor_TureJointIDA_rowGene_Rcof_Spearman","Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero","Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueJointIDA_EstJointIDA ","RMSE_TrueJointIDA_ImputeGene_Rcof","RMSE_TrueJointIDA_rowGene_Rcof","proportion_GeneToGene_Net","Cor_rowGene_ImputeGene_Rcof","shd_ture_estGrapph","shd_ture_randomGraph")
  list_related_variable = list()
  list_related_variable[["related_variable"]] = related_variable
  list_related_variable[["related_variable_colname"]] = related_variable_colname
  
  return(list_related_variable)
}


#***************************************************************************************************
# Aim of this function is similar as function 'Get_related_variable' with more summarized Statistical indicators
# This is an update version of 'Get_related_variable' function
#***************************************************************************************************
# @ - dt_IDA: the results of estimated IDA and jointIDA based on our causal graph
# @ - adj: adjacency matrix of true graph
# @ - g1: the simulated true causal graph
# @ - g2: the estimated causal graph based on our framework
# @ - g3: a random causal graph which is used for comparison
# @ - gene_index: index of genes in the simulated causal network
#***************************************************************************************************
Get_related_variable2<-function(dt_IDA,adj,g1,g2,g3,gene_index)
{
  
  #save correlation of IDA
  Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA,  use='complete.obs' ) 
  Cor_TureIDA_EstIDA_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA, use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson	= cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs')
  Cor_TureIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs', method="spearman")
  Cor_TureIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof, method="spearman")
  
  
  #save correlation of nonzero IDA
  ind_IDA =  which(abs(dt_IDA$true_IDA)>1e-10)
  Cor_TureIDA_EstIDA_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA[ind_IDA],  use='complete.obs' )
  Cor_TureIDA_EstIDA_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA[ind_IDA], use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs' )
  Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs', method="spearman" )
  Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA])
  Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA], method="spearman")
  
  #save RMSE
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA)
  RMSE_TrueIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueIDA_rowGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  
  #____________________________
  # correlation of joint IDA
  #____________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  # save correlation of joint IDA
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  Cor_TureJointIDA_EstJointIDA_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA, method="spearman")
  
  # modified by fyn 14/2/2020
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof,  method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureJointIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof,  method="spearman")
  
  
  #________________________________________________
  # correlation of joint IDA which are non-zero
  #________________________________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  ind_nonzero =  which(abs(as.numeric(dt_IDA$true_jointIDA))>1e-10)
  
  Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero])
  Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero], method="spearman")
  
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero])
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero],method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero])	
  Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero],  method="spearman")
  
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  RMSE_TrueJointIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueJointIDA_rowGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  # proportion of gene-gene connection
  gene_mat = adj[gene_index,gene_index]
  sum(gene_mat!=0)/nrow(gene_mat)^2
  
  proportion_GeneToGene_Net = sum(gene_mat!=0)/nrow(gene_mat)^2
  
  # correlation of regression coef derived from real expr vs imputed expr
  Cor_rowGene_ImputeGene_Rcof = cor(dt_IDA$Gene_rowdata_regression_Cof, dt_IDA$Gene_imputedata_regression_Cof)
  
  # shd_ture_estGrapph = shd(g1,g2)
  # shd_ture_randomGraph = shd(g1,g3)
  shd_ture_estGrapph = 0
  shd_ture_randomGraph = 0
  
  # correlation and MSE
  related_variable = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Spearman,Cor_TureIDA_ImputeGene_Rcof_Pearson,Cor_TureIDA_ImputeGene_Rcof_Spearman,	Cor_TureIDA_rowGene_Rcof_Pearson,Cor_TureIDA_rowGene_Rcof_Spearman,Cor_TureIDA_EstIDA_Pearson_Nonzero,	Cor_TureIDA_EstIDA_Spearman_Nonzero,	Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero,Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueIDA_EstIDA,	RMSE_TrueIDA_ImputeGene_Rcof,	RMSE_TrueIDA_rowGene_Rcof,	Cor_TureJointIDA_EstJointIDA_Pearson,	Cor_TureJointIDA_EstJointIDA_Spearman,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman,	Cor_TureJointIDA_rowGene_Rcof_Pearson,	Cor_TureJointIDA_rowGene_Rcof_Spearman,	Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero,	Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueJointIDA_EstJointIDA ,	RMSE_TrueJointIDA_ImputeGene_Rcof ,	RMSE_TrueJointIDA_rowGene_Rcof ,	proportion_GeneToGene_Net,	Cor_rowGene_ImputeGene_Rcof,	shd_ture_estGrapph,	shd_ture_randomGraph)
  related_variable_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Spearman","Cor_TureIDA_ImputeGene_Rcof_Pearson","Cor_TureIDA_ImputeGene_Rcof_Spearman","Cor_TureIDA_rowGene_Rcof_Pearson","Cor_TureIDA_rowGene_Rcof_Spearman","Cor_TureIDA_EstIDA_Pearson_Nonzero","Cor_TureIDA_EstIDA_Spearman_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_ImputeGene_Rcof","RMSE_TrueIDA_rowGene_Rcof","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Spearman","Cor_TureJointIDA_ImputeGene_Rcof_Pearson","Cor_TureJointIDA_ImputeGene_Rcof_Spearman","Cor_TureJointIDA_rowGene_Rcof_Pearson","Cor_TureJointIDA_rowGene_Rcof_Spearman","Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero","Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueJointIDA_EstJointIDA ","RMSE_TrueJointIDA_ImputeGene_Rcof","RMSE_TrueJointIDA_rowGene_Rcof","proportion_GeneToGene_Net","Cor_rowGene_ImputeGene_Rcof","shd_ture_estGrapph","shd_ture_randomGraph")
  list_related_variable = list()
  list_related_variable[["related_variable"]] = related_variable
  list_related_variable[["related_variable_colname"]] = related_variable_colname
  
  return(list_related_variable)
}



#***************************************************************************************************
# Aim of this function is to summarize the results of Lasso, MVR etc.
# Caculate the related variable,including correlation and RMSE when comparing other methods.(Lasso, MVR etc)
#***************************************************************************************************
# @ - dt_info_trueIDA: the data table includes IDA and jointIDA
# @ - est_coef_mat: results of other methods, including Lasso, MVR etc.
#***************************************************************************************************
Get_related_variable_compareOtherMethod<-function(dt_info_trueIDA,est_coef_mat)
{
  
  cor_trueIDA_MVLR_coef_Pearson = cor(dt_info_trueIDA$true_IDA, est_coef_mat$MVLR_coef,  use='complete.obs' ) 
  cor_trueIDA_Lasso_coef_Pearson = cor(dt_info_trueIDA$true_IDA, est_coef_mat$Lasso_coef,  use='complete.obs' ) 
  cor_trueIDA_Elasticnet_coef_Pearson = cor(dt_info_trueIDA$true_IDA, est_coef_mat$Elasticnet_coef,  use='complete.obs' ) 
  
  RMSE_trueIDA_MVLR_coef_Pearson = RMSE(dt_info_trueIDA$true_IDA, est_coef_mat$MVLR_coef) 
  RMSE_trueIDA_Lasso_coef_Pearson = RMSE(dt_info_trueIDA$true_IDA, est_coef_mat$Lasso_coef) 
  RMSE_trueIDA_Elasticnet_coef_Pearson = RMSE(dt_info_trueIDA$true_IDA, est_coef_mat$Elasticnet_coef) 
  
  # correlation and MSE
  related_variable = c(cor_trueIDA_MVLR_coef_Pearson, cor_trueIDA_Lasso_coef_Pearson, cor_trueIDA_Elasticnet_coef_Pearson, RMSE_trueIDA_MVLR_coef_Pearson, RMSE_trueIDA_Lasso_coef_Pearson, RMSE_trueIDA_Elasticnet_coef_Pearson)
  related_variable_colname = c("cor_trueIDA_MVLR_coef_Pearson", "cor_trueIDA_Lasso_coef_Pearson", 
                               "cor_trueIDA_Elasticnet_coef_Pearson", "RMSE_trueIDA_MVLR_coef_Pearson", 
                               "RMSE_trueIDA_Lasso_coef_Pearson", "RMSE_trueIDA_Elasticnet_coef_Pearson")
  
  list_related_variable = list()
  list_related_variable[["related_variable"]] = related_variable
  list_related_variable[["related_variable_colname"]] = related_variable_colname
  
  return(list_related_variable)
}



#***************************************************************************************************
# Aim of this function is to combine the results of Lasso, MVR etc.
#***************************************************************************************************
# @ - est_coef_eval_final: results of other methods, including Lasso, MVR etc.
#***************************************************************************************************
ReOrganize_MVLR_Lasso_EN_Res<-function(est_coef_eval_final)
{
  est_coef_eval_final = as.data.table(est_coef_eval_final)
  MVLR_res = est_coef_eval_final[est_coef_eval_final$Method=="MVLR",]
  MVLR_se = as.numeric(MVLR_res[,"rval_se"]$rval_se)
  MVLR_sp = as.numeric(MVLR_res[,"rval_sp"]$rval_sp)
  MVLR_ppv = as.numeric(MVLR_res[,"rval_ppv"]$rval_ppv)
  MVLR_npv = as.numeric(MVLR_res[,"rval_npv"]$rval_npv)
  MVLR_F1_score = as.numeric(2 * (MVLR_ppv * MVLR_se) / (MVLR_ppv + MVLR_se)) # Caculate the F1 score based on the precision and recall 
  MVLR_F2_score = as.numeric((1+2*2)*(MVLR_ppv * MVLR_se)/(2*2*MVLR_ppv+MVLR_se))
  MVLR_F0.5_score = as.numeric((1+0.5*0.5)*(MVLR_ppv * MVLR_se)/(0.5*0.5*MVLR_ppv+MVLR_se))
  
  Lasso_res = est_coef_eval_final[est_coef_eval_final$Method=="Lasso",]
  Lasso_se = as.numeric(Lasso_res[,"rval_se"]$rval_se)
  Lasso_sp = as.numeric(Lasso_res[,"rval_sp"]$rval_sp)
  Lasso_ppv = as.numeric(Lasso_res[,"rval_ppv"]$rval_ppv)
  Lasso_npv = as.numeric(Lasso_res[,"rval_npv"]$rval_npv)
  Lasso_F1_score = as.numeric(2 * (Lasso_ppv * Lasso_se) / (Lasso_ppv + Lasso_se)) # Caculate the F1 score based on the precision and recall 
  Lasso_F2_score = as.numeric((1+2*2)*(Lasso_ppv * Lasso_se)/(2*2*Lasso_ppv+Lasso_se))
  Lasso_F0.5_score = as.numeric((1+0.5*0.5)*(Lasso_ppv * Lasso_se)/(0.5*0.5*Lasso_ppv+Lasso_se))
  
  Elasticnet_res = est_coef_eval_final[est_coef_eval_final$Method=="Elasticnet",]
  Elasticnet_se = as.numeric(Elasticnet_res[,"rval_se"]$rval_se)
  Elasticnet_sp = as.numeric(Elasticnet_res[,"rval_sp"]$rval_sp)
  Elasticnet_ppv = as.numeric(Elasticnet_res[,"rval_ppv"]$rval_ppv)
  Elasticnet_npv = as.numeric(Elasticnet_res[,"rval_npv"]$rval_npv)
  Elasticnet_F1_score = as.numeric(2 * (Elasticnet_ppv * Elasticnet_se) / (Elasticnet_ppv + Elasticnet_se)) # Caculate the F1 score based on the precision and recall 
  Elasticnet_F2_score = as.numeric((1+2*2)*(Elasticnet_ppv * Elasticnet_se)/(2*2*Elasticnet_ppv+Elasticnet_se))
  Elasticnet_F0.5_score = as.numeric((1+0.5*0.5)*(Elasticnet_ppv * Elasticnet_se)/(0.5*0.5*Elasticnet_ppv+Elasticnet_se))
  
  
  
  
  related_variable = c(MVLR_se,MVLR_sp,MVLR_ppv,MVLR_npv,MVLR_F1_score,MVLR_F2_score,MVLR_F0.5_score,
                       Lasso_se,Lasso_sp,Lasso_ppv,Lasso_npv,Lasso_F1_score,Lasso_F2_score,Lasso_F0.5_score,
                       Elasticnet_se,Elasticnet_sp,Elasticnet_ppv,Elasticnet_npv,Elasticnet_F1_score,Elasticnet_F2_score,Elasticnet_F0.5_score)
  
  related_variable_colname = c("MVLR_se","MVLR_sp","MVLR_ppv","MVLR_npv","MVLR_F1_score","MVLR_F2_score","MVLR_F0.5_score",
                               "Lasso_se","Lasso_sp","Lasso_ppv","Lasso_npv","Lasso_F1_score","Lasso_F2_score","Lasso_F0.5_score",
                               "Elasticnet_se","Elasticnet_sp","Elasticnet_ppv","Elasticnet_npv","Elasticnet_F1_score","Elasticnet_F2_score","Elasticnet_F0.5_score")
  
  list_MVLR_Lasso_EN_Res = list()
  list_MVLR_Lasso_EN_Res[["related_variable"]] = related_variable
  list_MVLR_Lasso_EN_Res[["related_variable_colname"]] = related_variable_colname
  
  return(list_MVLR_Lasso_EN_Res)
  
}



#***************************************************************************************************
# Aim of this function is to estimate the performance of detecting mediating genes
#***************************************************************************************************
# @ - ifmediator_gene_true: results of mediating genes based on true causal graph
# @ - ifmediator_gene_est: results of mediating genes based on estimated causal graph
# @ - genes_num: the number of genes
#***************************************************************************************************
Performance_mediatorDetect_FromGraphStructure<-function(ifmediator_gene_true,ifmediator_gene_est,genes_num)
{
  ######################################################################################
  # Find the mediator in the gene-outcome graph by using graph structure.
  ifmediator_gene_true_std = Standardize_Mediator_result(ifmediator_gene_true)
  ifmediator_gene_est_std = Standardize_Mediator_result(ifmediator_gene_est)
  
  ifmediator_gene_true = ifmediator_gene_true[1:genes_num]
  ifmediator_gene_est = ifmediator_gene_est[1:genes_num]
  
  tab <- table(ifmediator_gene_est_std,ifmediator_gene_true_std)[2:1,2:1]
  mediator_rval <- epi.tests(tab, conf.level = 0.95)
  mediator_rval_aprev = mediator_rval$detail[1,2] # apparent prevalence.
  mediator_rval_tprev = mediator_rval$detail[2,2]  # true prevalence.
  mediator_rval_se = mediator_rval$detail[3,2]   # Sensitivity
  mediator_rval_sp = mediator_rval$detail[4,2]   # Specificity
  mediator_rval_ppv = mediator_rval$detail[9,2]   # Positive predictive value
  mediator_rval_npv = mediator_rval$detail[10,2]   # Negative predictive value
  mediator_rval_plr = mediator_rval$detail[11,2]   # Positive likelihood ratio
  mediator_rval_nlr = mediator_rval$detail[12,2]   # Negative likelihood ratio
  
  mediator_rval_combine = c(mediator_rval_aprev,mediator_rval_tprev,mediator_rval_se,mediator_rval_sp,mediator_rval_ppv,mediator_rval_npv,mediator_rval_plr,mediator_rval_nlr)
  mediator_rval_combine_colname = c("mediator_rval_aprev","mediator_rval_tprev","mediator_rval_se","mediator_rval_sp","mediator_rval_ppv","mediator_rval_npv","mediator_rval_plr","mediator_rval_nlr")
  
  list_mediator_rval_combine = list()
  list_mediator_rval_combine[["mediator_rval_combine"]] = mediator_rval_combine
  list_mediator_rval_combine[["mediator_rval_combine_colname"]] = mediator_rval_combine_colname
  
  return(list_mediator_rval_combine)
}



#***************************************************************************************************
# Aim of this function is to estimate the performance of detecting causal genes,including both direct and indirect genes
#***************************************************************************************************
# @ - gene_cause_true : results of causal genes based on true causal graph
# @ - gene_cause_est: results of causal genes based on estimated causal graph
#***************************************************************************************************
Performance_CausalGene<-function(gene_cause_est,gene_cause_true)
{
  tab <- table(gene_cause_est,gene_cause_true)[2:1,2:1]
  causalgene_rval <- epi.tests(tab, conf.level = 0.95)
  causalgene_rval_aprev = causalgene_rval$detail[1,2] # apparent prevalence.
  causalgene_rval_tprev = causalgene_rval$detail[2,2]  # true prevalence.
  causalgene_rval_se = causalgene_rval$detail[3,2]   # Sensitivity
  causalgene_rval_sp = causalgene_rval$detail[4,2]   # Specificity
  causalgene_rval_ppv = causalgene_rval$detail[9,2]   # Positive predictive value
  causalgene_rval_npv = causalgene_rval$detail[10,2]   # Negative predictive value
  causalgene_rval_plr = causalgene_rval$detail[11,2]   # Positive likelihood ratio
  causalgene_rval_nlr = causalgene_rval$detail[12,2]   # Negative likelihood ratio
  
  causalgene_rval_combine = c(causalgene_rval_aprev,causalgene_rval_tprev,causalgene_rval_se,causalgene_rval_sp,causalgene_rval_ppv,causalgene_rval_npv,causalgene_rval_plr,causalgene_rval_nlr)
  causalgene_rval_combine_colname = c("causalgene_rval_aprev","causalgene_rval_tprev","causalgene_rval_se","causalgene_rval_sp","causalgene_rval_ppv","causalgene_rval_npv","causalgene_rval_plr","causalgene_rval_nlr")
  
  list_causalgene_rval_combine = list()
  list_causalgene_rval_combine[["causalgene_rval_combine"]] = causalgene_rval_combine
  list_causalgene_rval_combine[["causalgene_rval_combine_colname"]] = causalgene_rval_combine_colname
  
  return(list_causalgene_rval_combine)
}



#***************************************************************************************************
# Aim of this function is to estimate the performance of learning gene-pheno causal network
#***************************************************************************************************
# @ - adj: adjacency matrix of true causal graph
# @ - gene_index: index of genes in the simulated causal network
# @ - n3: the total number of nodes in the causal network
# @ - pcSimple.fi: the results of gene-pheno causal network
#***************************************************************************************************
PC_Simple_Performance<-function(adj,gene_index,n3,pcSimple.fi)
{
  #_________________________________________________
  # this part examines the performance of PC simple in detecting direct causal variants
  #_____________________________________________________
  Pheno_vec = adj[gene_index,n3]
  
  gene_list = attributes(pcSimple.fit$G)
  pcsimple_index = attributes (which(pcSimple.fit$G==TRUE))$names # fyn: gene index with true relations to phenotype
  pcsimple_index = as.numeric(pcsimple_index) 
  real_directCausal = gene_index[Pheno_vec!=0]
  
  cell_11 = length(intersect(pcsimple_index, real_directCausal))
  cell_12 = length(pcsimple_index) - cell_11
  cell_21 = length(real_directCausal) - cell_11
  cell_22 = MyData[x,2] - cell_11 - cell_12 - cell_21  # MyData[x,2] is gene_num
  
  PCsimple_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  
  
  return(PCsimple_mat)
}


#***************************************************************************************************
# updated version of 'PC_Simple_Performance' function
#***************************************************************************************************
# @ - adj: adjacency matrix of true causal graph
# @ - gene_index: index of genes in the simulated causal network
# @ - n3: the total number of nodes in the causal network
# @ - pcSimple.fi: the results of gene-pheno causal network
#***************************************************************************************************
PC_Simple_Performance2<-function(adj,gene_index,n3,pcSimple.fi,genes_num)
{
  #_________________________________________________
  # this part examines the performance of PC simple in detecting direct causal variants
  #_____________________________________________________
  Pheno_vec = adj[gene_index,n3]
  
  gene_list = attributes(pcSimple.fit$G)
  pcsimple_index = attributes (which(pcSimple.fit$G==TRUE))$names # fyn: gene index with true relations to phenotype
  pcsimple_index = as.numeric(pcsimple_index) 
  real_directCausal = gene_index[Pheno_vec!=0]
  
  cell_11 = length(intersect(pcsimple_index, real_directCausal))
  cell_12 = length(pcsimple_index) - cell_11
  cell_21 = length(real_directCausal) - cell_11
  cell_22 = genes_num - cell_11 - cell_12 - cell_21  # MyData[x,2] is gene_num
  
  PCsimple_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  PCsimple_rval <- epi.tests(PCsimple_mat, conf.level = 0.95)
  PCsimple_rval_aprev = PCsimple_rval$detail[1,2] # apparent prevalence.
  PCsimple_rval_tprev = PCsimple_rval$detail[2,2]  # true prevalence.
  PCsimple_rval_se = PCsimple_rval$detail[3,2]   # Sensitivity
  PCsimple_rval_sp = PCsimple_rval$detail[4,2]   # Specificity
  PCsimple_rval_ppv = PCsimple_rval$detail[9,2]   # Positive predictive value
  PCsimple_rval_npv = PCsimple_rval$detail[10,2]   # Negative predictive value
  PCsimple_rval_plr = PCsimple_rval$detail[11,2]   # Positive likelihood ratio
  PCsimple_rval_nlr = PCsimple_rval$detail[12,2]   # Negative likelihood ratio
  
  PCsimple_rval_combine = c(PCsimple_rval_aprev,PCsimple_rval_tprev,PCsimple_rval_se,PCsimple_rval_sp,PCsimple_rval_ppv,PCsimple_rval_npv,PCsimple_rval_plr,PCsimple_rval_nlr)
  PCsimple_rval_combine_colname = c("PCsimple_rval_aprev","PCsimple_rval_tprev","PCsimple_rval_se","PCsimple_rval_sp","PCsimple_rval_ppv","PCsimple_rval_npv","PCsimple_rval_plr","PCsimple_rval_nlr")
  
  list_PCsimple_rval_combine = list()
  list_PCsimple_rval_combine[["PCsimple_rval_combine"]] = PCsimple_rval_combine
  list_PCsimple_rval_combine[["PCsimple_rval_combine_colname"]] = PCsimple_rval_combine_colname
  
  
  return(list_PCsimple_rval_combine)
}


#***************************************************************************************************
# Aim of this function is to estimate the performance of univariate test based on the imputed gene expression data
#***************************************************************************************************
# @ - dt_IDA: data table includes the results of IDA and jointIDA
# @ - MyData: data includes the information of scenarios
# @ - Pheno_vec: information includes direct causal genes with phenotype
#***************************************************************************************************
Impute_Performance<-function(dt_IDA,MyData,Pheno_vec)
{
  # this part examines the performance of PrediXcan (TWAS) in detecting direct causal variants
  #_____________________________________________________
  impExpr_ind = which( dt_IDA$Gene_imputedata_regression_pval <= 0.001)
  real_ind = which(Pheno_vec!=0) 
  cell_11 = length(  intersect(impExpr_ind, real_ind)  )
  cell_12 = length(impExpr_ind) - cell_11
  cell_21 = length(real_ind) - cell_11
 
  cell_22 = MyData[x,2] - cell_11 - cell_12 - cell_21
  ImpExpr_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  return(ImpExpr_mat)
}



#***************************************************************************************************
# This is an updated version of function 'Impute_Performance'
# The main difference between this function and the previous one is we give the gene_num to the function directly
#***************************************************************************************************
# @ - dt_IDA: data table includes the results of IDA and jointIDA
# @ - gene_num: number of genes
# @ - Pheno_vec: information includes direct causal genes with phenotype
#***************************************************************************************************
Impute_Performance2<-function(dt_IDA,gene_num,Pheno_vec)
{
  # this part examines the performance of PrediXcan (TWAS) in detecting direct causal variants
  #_____________________________________________________
  impExpr_ind = which( dt_IDA$Gene_imputedata_regression_pval <= 0.001)
  real_ind = which(Pheno_vec!=0) 
  cell_11 = length(  intersect(impExpr_ind, real_ind)  )
  cell_12 = length(impExpr_ind) - cell_11
  cell_21 = length(real_ind ) - cell_11
  
  cell_22 = gene_num - cell_11 - cell_12 - cell_21
  ImpExpr_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  return(ImpExpr_mat)
}



#***************************************************************************************************
# This is an updated version of function 'Impute_Performance2'
# add the epi.test for the matrix
#***************************************************************************************************
# @ - dt_IDA: data table includes the results of IDA and jointIDA
# @ - gene_num: number of genes
# @ - Pheno_vec: information includes direct causal genes with phenotype
#***************************************************************************************************
Impute_Performance3<-function(dt_IDA,gene_num,Pheno_vec)
{
  # this part examines the performance of PrediXcan (TWAS) in detecting direct causal variants
  #_____________________________________________________
  impExpr_ind = which( dt_IDA$Gene_imputedata_regression_pval <= 0.001)
  real_ind = which(Pheno_vec!=0) 
  cell_11 = length(  intersect(impExpr_ind, real_ind)  )
  cell_12 = length(impExpr_ind) - cell_11
  cell_21 = length(real_ind ) - cell_11
  
  cell_22 = gene_num - cell_11 - cell_12 - cell_21
  ImpExpr_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  ImpExpr_rval <- epi.tests(ImpExpr_mat, conf.level = 0.95)
  ImpExpr_rval_aprev = ImpExpr_rval$detail[1,2] # apparent prevalence.
  ImpExpr_rval_tprev = ImpExpr_rval$detail[2,2]  # true prevalence.
  ImpExpr_rval_se = ImpExpr_rval$detail[3,2]   # Sensitivity
  ImpExpr_rval_sp = ImpExpr_rval$detail[4,2]   # Specificity
  ImpExpr_rval_ppv = ImpExpr_rval$detail[9,2]   # Positive predictive value
  ImpExpr_rval_npv = ImpExpr_rval$detail[10,2]   # Negative predictive value
  ImpExpr_rval_plr = ImpExpr_rval$detail[11,2]   # Positive likelihood ratio
  ImpExpr_rval_nlr = ImpExpr_rval$detail[12,2]   # Negative likelihood ratio
  
  ImpExpr_rval_combine = c(ImpExpr_rval_aprev,ImpExpr_rval_tprev,ImpExpr_rval_se,ImpExpr_rval_sp,ImpExpr_rval_ppv,ImpExpr_rval_npv,ImpExpr_rval_plr,ImpExpr_rval_nlr)
  ImpExpr_rval_combine_colname = c("ImpExpr_rval_aprev","ImpExpr_rval_tprev","ImpExpr_rval_se","ImpExpr_rval_sp","ImpExpr_rval_ppv","ImpExpr_rval_npv","ImpExpr_rval_plr","ImpExpr_rval_nlr")
  
  list_ImpExpr_rval_combine = list()
  list_ImpExpr_rval_combine[["ImpExpr_rval_combine"]] = ImpExpr_rval_combine
  list_ImpExpr_rval_combine[["ImpExpr_rval_combine_colname"]] = ImpExpr_rval_combine_colname
  
  return(list_ImpExpr_rval_combine)
}


#***************************************************************************************************
# This is an updated version of function 'Impute_Performance3'
### Added by fyn. 2023.10.4. Using the results from univariate test directly to simplify the code
#***************************************************************************************************
Impute_Performance4<-function(univariate_test_Matrix,gene_num,Pheno_vec)
{
  # this part examines the performance of PrediXcan (TWAS) in detecting direct causal variants
  #_____________________________________________________
  impExpr_ind = which( univariate_test_Matrix$impgene_pval <= 0.001)
  real_ind = which(Pheno_vec!=0) 
  cell_11 = length(  intersect(impExpr_ind, real_ind)  )
  cell_12 = length(impExpr_ind) - cell_11
  cell_21 = length(real_ind) - cell_11
  
  cell_22 = MyData[x,2] - cell_11 - cell_12 - cell_21
  ImpExpr_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  return(ImpExpr_mat)
}



Get_epi_fisher_Summary<-function(PCsimple_mat,ImpExpr_mat)
{
  epi.tests_PCsimple = epi.tests(PCsimple_mat)
  fisher.test_PCsimple = fisher.test(PCsimple_mat)
  
  epi.tests_ImpExpr = epi.tests(ImpExpr_mat)
  fisher.test_ImpExpr = fisher.test(ImpExpr_mat)
  
  epi_PCsimple_aprev = epi.tests_PCsimple$detail[1,2] # apparent prevalence.
  epi_PCsimple_tprev = epi.tests_PCsimple$detail[2,2]  # true prevalence.
  epi_PCsimple_se = epi.tests_PCsimple$detail[3,2]   # Sensitivity
  epi_PCsimple_sp = epi.tests_PCsimple$detail[4,2]   # Specificity
  epi_PCsimple_ppv = epi.tests_PCsimple$detail[9,2]   # Positive predictive value
  epi_PCsimple_npv = epi.tests_PCsimple$detail[10,2]   # Negative predictive value
  epi_PCsimple_plr = epi.tests_PCsimple$detail[11,2]   # Positive likelihood ratio
  epi_PCsimple_nlr = epi.tests_PCsimple$detail[12,2]   # Negative likelihood ratio
  
  epi_ImpExpr_aprev = epi.tests_ImpExpr$detail[1,2] # apparent prevalence.
  epi_ImpExpr_tprev = epi.tests_ImpExpr$detail[2,2]  # true prevalence.
  epi_ImpExpr_se = epi.tests_ImpExpr$detail[3,2]   # Sensitivity
  epi_ImpExpr_sp = epi.tests_ImpExpr$detail[4,2]   # Specificity
  epi_ImpExpr_ppv = epi.tests_ImpExpr$detail[9,2]   # Positive predictive value
  epi_ImpExpr_npv = epi.tests_ImpExpr$detail[10,2]   # Negative predictive value
  epi_ImpExpr_plr = epi.tests_ImpExpr$detail[11,2]   # Positive likelihood ratio
  epi_ImpExpr_nlr = epi.tests_ImpExpr$detail[12,2]   # Negative likelihood ratio
  
  # ##extract each data in epi.test
  # epi_PCsimple_aprev_est = epi_PCsimple_aprev[1,1]
  # epi_PCsimple_aprev_lower = epi_PCsimple_aprev[1,2]
  # epi_PCsimple_aprev_upper = epi_PCsimple_aprev[1,3]
  # epi_PCsimple_tprev_est = epi_PCsimple_tprev[1,1]
  # epi_PCsimple_tprev_lower = epi_PCsimple_tprev[1,2]
  # epi_PCsimple_tprev_upper = epi_PCsimple_tprev[1,3]
  # epi_PCsimple_se_est = epi_PCsimple_se[1,1]
  # epi_PCsimple_se_lower = epi_PCsimple_se[1,2]
  # epi_PCsimple_se_upper = epi_PCsimple_se[1,3]
  # epi_PCsimple_sp_est = epi_PCsimple_sp[1,1]
  # epi_PCsimple_sp_lower = epi_PCsimple_sp[1,2]
  # epi_PCsimple_sp_upper = epi_PCsimple_sp[1,3]
  # epi_PCsimple_ppv_est = epi_PCsimple_ppv[1,1]
  # epi_PCsimple_ppv_lower = epi_PCsimple_ppv[1,2]
  # epi_PCsimple_ppv_upper = epi_PCsimple_ppv[1,3]
  # epi_PCsimple_npv_est = epi_PCsimple_npv[1,1]
  # epi_PCsimple_npv_lower = epi_PCsimple_npv[1,2]
  # epi_PCsimple_npv_upper = epi_PCsimple_npv[1,3]
  # epi_PCsimple_plr_est = epi_PCsimple_plr[1,1]
  # epi_PCsimple_plr_lower = epi_PCsimple_plr[1,2]
  # epi_PCsimple_plr_upper = epi_PCsimple_plr[1,3]
  # epi_PCsimple_nlr_est = epi_PCsimple_nlr[1,1]
  # epi_PCsimple_nlr_lower = epi_PCsimple_nlr[1,2]
  # epi_PCsimple_nlr_upper = epi_PCsimple_nlr[1,3]
  # 
  # epi_ImpExpr_aprev_est = epi_ImpExpr_aprev[1,1]
  # epi_ImpExpr_aprev_lower = epi_ImpExpr_aprev[1,2]
  # epi_ImpExpr_aprev_upper = epi_ImpExpr_aprev[1,3]
  # epi_ImpExpr_tprev_est = epi_ImpExpr_tprev[1,1]
  # epi_ImpExpr_tprev_lower = epi_ImpExpr_tprev[1,2]
  # epi_ImpExpr_tprev_upper = epi_ImpExpr_tprev[1,3]
  # epi_ImpExpr_se_est = epi_ImpExpr_se[1,1]
  # epi_ImpExpr_se_lower = epi_ImpExpr_se[1,2]
  # epi_ImpExpr_se_upper = epi_ImpExpr_se[1,3]
  # epi_ImpExpr_sp_est = epi_ImpExpr_sp[1,1]
  # epi_ImpExpr_sp_lower = epi_ImpExpr_sp[1,2]
  # epi_ImpExpr_sp_upper = epi_ImpExpr_sp[1,3]
  # epi_ImpExpr_ppv_est = epi_ImpExpr_ppv[1,1]
  # epi_ImpExpr_ppv_lower = epi_ImpExpr_ppv[1,2]
  # epi_ImpExpr_ppv_upper = epi_ImpExpr_ppv[1,3]
  # epi_ImpExpr_npv_est = epi_ImpExpr_npv[1,1]
  # epi_ImpExpr_npv_lower = epi_ImpExpr_npv[1,2]
  # epi_ImpExpr_npv_upper = epi_ImpExpr_npv[1,3]
  # epi_ImpExpr_plr_est = epi_ImpExpr_plr[1,1]
  # epi_ImpExpr_plr_lower = epi_ImpExpr_plr[1,2]
  # epi_ImpExpr_plr_upper = epi_ImpExpr_plr[1,3]
  # epi_ImpExpr_nlr_est = epi_ImpExpr_nlr[1,1]
  # epi_ImpExpr_nlr_lower = epi_ImpExpr_nlr[1,2]
  # epi_ImpExpr_nlr_upper = epi_ImpExpr_nlr[1,3]
  # 
  # 
  # fisher_PCsimple_pvalue = fisher.test(PCsimple_mat)$p.value
  # fisher_PCsimple_95PCI_1 = fisher.test(PCsimple_mat)$conf.int[1] # 95 percent confidence interval
  # fisher_PCsimple_95PCI_2 = fisher.test(PCsimple_mat)$conf.int[2] # 95 percent confidence interval
  # fisher_PCsimple_OR = as.numeric(fisher.test(PCsimple_mat)$estimate) # odds ratio
  # 
  # fisher_ImpExpr_pvalue = fisher.test(ImpExpr_mat)$p.value
  # fisher_ImpExpr_95PCI_1 = fisher.test(ImpExpr_mat)$conf.int[1] # 95 percent confidence interval
  # fisher_ImpExpr_95PCI_2 = fisher.test(ImpExpr_mat)$conf.int[2] # 95 percent confidence interval
  # fisher_ImpExpr_OR = as.numeric(fisher.test(ImpExpr_mat)$estimate) # odds ratio
  
  epi_fisher_testresult_colname = c("epi_PCsimple_aprev","epi_PCsimple_tprev","epi_PCsimple_se",
                                    "epi_PCsimple_sp","epi_PCsimple_ppv","epi_PCsimple_npv",
                                    "epi_PCsimple_plr","epi_PCsimple_nlr","epi_ImpExpr_aprev",
                                    "epi_ImpExpr_tprev","epi_ImpExpr_se","epi_ImpExpr_sp",
                                    "epi_ImpExpr_ppv","epi_ImpExpr_npv","epi_ImpExpr_plr","epi_ImpExpr_nlr")
  epi_fisher_testresult = c(epi_PCsimple_aprev,epi_PCsimple_tprev,epi_PCsimple_se,
                            epi_PCsimple_sp,epi_PCsimple_ppv,epi_PCsimple_npv,
                            epi_PCsimple_plr,epi_PCsimple_nlr,epi_ImpExpr_aprev,
                            epi_ImpExpr_tprev,epi_ImpExpr_se,epi_ImpExpr_sp,epi_ImpExpr_ppv,
                            epi_ImpExpr_npv,epi_ImpExpr_plr,epi_ImpExpr_nlr)
  
  
  list_epi_fisher_testresult = list()
  list_epi_fisher_testresult[["epi_fisher_testresult_colname"]] = epi_fisher_testresult_colname
  list_epi_fisher_testresult[["epi_fisher_testresult"]] = epi_fisher_testresult
  
  return(list_epi_fisher_testresult)
}

Get_epi_fisher_Summary_pre<-function(PCsimple_mat,ImpExpr_mat)
{
  epi.tests_PCsimple = epi.tests(PCsimple_mat)
  fisher.test_PCsimple = fisher.test(PCsimple_mat)
  
  epi.tests_ImpExpr = epi.tests(ImpExpr_mat)
  fisher.test_ImpExpr = fisher.test(ImpExpr_mat)
  
  epi_PCsimple_aprev = epi.tests_PCsimple$rval['aprev'][[1]] # apparent prevalence.
  epi_PCsimple_tprev = epi.tests_PCsimple$rval['tprev'][[1]]  # true prevalence.
  epi_PCsimple_se = epi.tests_PCsimple$rval['se'][[1]]   # Sensitivity
  epi_PCsimple_sp = epi.tests_PCsimple$rval['sp'][[1]]   # Specificity
  epi_PCsimple_ppv = epi.tests_PCsimple$rval['ppv'][[1]]   # Positive predictive value
  epi_PCsimple_npv = epi.tests_PCsimple$rval['npv'][[1]]   # Negative predictive value
  epi_PCsimple_plr = epi.tests_PCsimple$rval['plr'][[1]]   # Positive likelihood ratio
  epi_PCsimple_nlr = epi.tests_PCsimple$rval['nlr'][[1]]   # Negative likelihood ratio
  
  epi_ImpExpr_aprev = epi.tests_ImpExpr$rval['aprev'][[1]] # apparent prevalence.
  epi_ImpExpr_tprev = epi.tests_ImpExpr$rval['tprev'][[1]]  # true prevalence.
  epi_ImpExpr_se = epi.tests_ImpExpr$rval['se'][[1]]   # Sensitivity
  epi_ImpExpr_sp = epi.tests_ImpExpr$rval['sp'][[1]]   # Specificity
  epi_ImpExpr_ppv = epi.tests_ImpExpr$rval['ppv'][[1]]   # Positive predictive value
  epi_ImpExpr_npv = epi.tests_ImpExpr$rval['npv'][[1]]   # Negative predictive value
  epi_ImpExpr_plr = epi.tests_ImpExpr$rval['plr'][[1]]   # Positive likelihood ratio
  epi_ImpExpr_nlr = epi.tests_ImpExpr$rval['nlr'][[1]]   # Negative likelihood ratio
  
  ##extract each data in epi.test
  epi_PCsimple_aprev_est = epi_PCsimple_aprev[1,1]
  epi_PCsimple_aprev_lower = epi_PCsimple_aprev[1,2]
  epi_PCsimple_aprev_upper = epi_PCsimple_aprev[1,3]
  epi_PCsimple_tprev_est = epi_PCsimple_tprev[1,1]
  epi_PCsimple_tprev_lower = epi_PCsimple_tprev[1,2]
  epi_PCsimple_tprev_upper = epi_PCsimple_tprev[1,3]
  epi_PCsimple_se_est = epi_PCsimple_se[1,1]
  epi_PCsimple_se_lower = epi_PCsimple_se[1,2]
  epi_PCsimple_se_upper = epi_PCsimple_se[1,3]
  epi_PCsimple_sp_est = epi_PCsimple_sp[1,1]
  epi_PCsimple_sp_lower = epi_PCsimple_sp[1,2]
  epi_PCsimple_sp_upper = epi_PCsimple_sp[1,3]
  epi_PCsimple_ppv_est = epi_PCsimple_ppv[1,1]
  epi_PCsimple_ppv_lower = epi_PCsimple_ppv[1,2]
  epi_PCsimple_ppv_upper = epi_PCsimple_ppv[1,3]
  epi_PCsimple_npv_est = epi_PCsimple_npv[1,1]
  epi_PCsimple_npv_lower = epi_PCsimple_npv[1,2]
  epi_PCsimple_npv_upper = epi_PCsimple_npv[1,3]
  epi_PCsimple_plr_est = epi_PCsimple_plr[1,1]
  epi_PCsimple_plr_lower = epi_PCsimple_plr[1,2]
  epi_PCsimple_plr_upper = epi_PCsimple_plr[1,3]
  epi_PCsimple_nlr_est = epi_PCsimple_nlr[1,1]
  epi_PCsimple_nlr_lower = epi_PCsimple_nlr[1,2]
  epi_PCsimple_nlr_upper = epi_PCsimple_nlr[1,3]
  
  epi_ImpExpr_aprev_est = epi_ImpExpr_aprev[1,1]
  epi_ImpExpr_aprev_lower = epi_ImpExpr_aprev[1,2]
  epi_ImpExpr_aprev_upper = epi_ImpExpr_aprev[1,3]
  epi_ImpExpr_tprev_est = epi_ImpExpr_tprev[1,1]
  epi_ImpExpr_tprev_lower = epi_ImpExpr_tprev[1,2]
  epi_ImpExpr_tprev_upper = epi_ImpExpr_tprev[1,3]
  epi_ImpExpr_se_est = epi_ImpExpr_se[1,1]
  epi_ImpExpr_se_lower = epi_ImpExpr_se[1,2]
  epi_ImpExpr_se_upper = epi_ImpExpr_se[1,3]
  epi_ImpExpr_sp_est = epi_ImpExpr_sp[1,1]
  epi_ImpExpr_sp_lower = epi_ImpExpr_sp[1,2]
  epi_ImpExpr_sp_upper = epi_ImpExpr_sp[1,3]
  epi_ImpExpr_ppv_est = epi_ImpExpr_ppv[1,1]
  epi_ImpExpr_ppv_lower = epi_ImpExpr_ppv[1,2]
  epi_ImpExpr_ppv_upper = epi_ImpExpr_ppv[1,3]
  epi_ImpExpr_npv_est = epi_ImpExpr_npv[1,1]
  epi_ImpExpr_npv_lower = epi_ImpExpr_npv[1,2]
  epi_ImpExpr_npv_upper = epi_ImpExpr_npv[1,3]
  epi_ImpExpr_plr_est = epi_ImpExpr_plr[1,1]
  epi_ImpExpr_plr_lower = epi_ImpExpr_plr[1,2]
  epi_ImpExpr_plr_upper = epi_ImpExpr_plr[1,3]
  epi_ImpExpr_nlr_est = epi_ImpExpr_nlr[1,1]
  epi_ImpExpr_nlr_lower = epi_ImpExpr_nlr[1,2]
  epi_ImpExpr_nlr_upper = epi_ImpExpr_nlr[1,3]
  
  
  fisher_PCsimple_pvalue = fisher.test(PCsimple_mat)$p.value
  fisher_PCsimple_95PCI_1 = fisher.test(PCsimple_mat)$conf.int[1] # 95 percent confidence interval
  fisher_PCsimple_95PCI_2 = fisher.test(PCsimple_mat)$conf.int[2] # 95 percent confidence interval
  fisher_PCsimple_OR = as.numeric(fisher.test(PCsimple_mat)$estimate) # odds ratio
  
  fisher_ImpExpr_pvalue = fisher.test(ImpExpr_mat)$p.value
  fisher_ImpExpr_95PCI_1 = fisher.test(ImpExpr_mat)$conf.int[1] # 95 percent confidence interval
  fisher_ImpExpr_95PCI_2 = fisher.test(ImpExpr_mat)$conf.int[2] # 95 percent confidence interval
  fisher_ImpExpr_OR = as.numeric(fisher.test(ImpExpr_mat)$estimate) # odds ratio
  
  epi_fisher_testresult_colname = c("epi_PCsimple_aprev_est ","epi_PCsimple_tprev_est ","epi_PCsimple_se_est ","epi_PCsimple_sp_est ","epi_PCsimple_ppv_est ","epi_PCsimple_npv_est ","epi_PCsimple_plr_est ","epi_PCsimple_nlr_est ","epi_ImpExpr_aprev_est ","epi_ImpExpr_tprev_est ","epi_ImpExpr_se_est ","epi_ImpExpr_sp_est ","epi_ImpExpr_ppv_est ","epi_ImpExpr_npv_est ","epi_ImpExpr_plr_est ","epi_ImpExpr_nlr_est","fisher_PCsimple_pvalue ","fisher_PCsimple_95PCI_1 ","fisher_PCsimple_95PCI_2 ","fisher_PCsimple_OR","fisher_ImpExpr_pvalue","fisher_ImpExpr_95PCI_1","fisher_ImpExpr_95PCI_2","fisher_ImpExpr_OR","epi_PCsimple_aprev_lower ","epi_PCsimple_aprev_upper","epi_PCsimple_tprev_lower ","epi_PCsimple_tprev_upper","epi_PCsimple_se_lower ","epi_PCsimple_se_upper ","epi_PCsimple_sp_lower ","epi_PCsimple_sp_upper ","epi_PCsimple_ppv_lower ","epi_PCsimple_ppv_upper ","epi_PCsimple_npv_lower ","epi_PCsimple_npv_upper ","epi_PCsimple_plr_lower ","epi_PCsimple_plr_upper","epi_PCsimple_nlr_lower ","epi_PCsimple_nlr_upper ","epi_ImpExpr_aprev_lower","epi_ImpExpr_aprev_upper","epi_ImpExpr_tprev_lower","epi_ImpExpr_tprev_upper","epi_ImpExpr_se_lower","epi_ImpExpr_se_upper","epi_ImpExpr_sp_lower","epi_ImpExpr_sp_upper","epi_ImpExpr_ppv_lower","epi_ImpExpr_ppv_upper","epi_ImpExpr_npv_lower","epi_ImpExpr_npv_upper","epi_ImpExpr_plr_lower","epi_ImpExpr_plr_upper","epi_ImpExpr_nlr_lower","epi_ImpExpr_nlr_upper")
  epi_fisher_testresult = c(epi_PCsimple_aprev_est ,	epi_PCsimple_tprev_est ,	epi_PCsimple_se_est ,	epi_PCsimple_sp_est ,	epi_PCsimple_ppv_est ,	epi_PCsimple_npv_est ,	epi_PCsimple_plr_est ,	epi_PCsimple_nlr_est ,	epi_ImpExpr_aprev_est ,	epi_ImpExpr_tprev_est ,	epi_ImpExpr_se_est ,	epi_ImpExpr_sp_est ,	epi_ImpExpr_ppv_est ,	epi_ImpExpr_npv_est ,	epi_ImpExpr_plr_est ,	epi_ImpExpr_nlr_est,	fisher_PCsimple_pvalue ,	fisher_PCsimple_95PCI_1 ,	fisher_PCsimple_95PCI_2 ,	fisher_PCsimple_OR,	fisher_ImpExpr_pvalue,	fisher_ImpExpr_95PCI_1,	fisher_ImpExpr_95PCI_2,	fisher_ImpExpr_OR,	epi_PCsimple_aprev_lower ,	epi_PCsimple_aprev_upper,	epi_PCsimple_tprev_lower ,	epi_PCsimple_tprev_upper,	epi_PCsimple_se_lower ,	epi_PCsimple_se_upper ,	epi_PCsimple_sp_lower ,	epi_PCsimple_sp_upper ,	epi_PCsimple_ppv_lower ,	epi_PCsimple_ppv_upper ,	epi_PCsimple_npv_lower ,	epi_PCsimple_npv_upper ,	epi_PCsimple_plr_lower ,	epi_PCsimple_plr_upper,	epi_PCsimple_nlr_lower ,	epi_PCsimple_nlr_upper ,	epi_ImpExpr_aprev_lower,	epi_ImpExpr_aprev_upper,	epi_ImpExpr_tprev_lower,	epi_ImpExpr_tprev_upper,	epi_ImpExpr_se_lower,	epi_ImpExpr_se_upper,	epi_ImpExpr_sp_lower,	epi_ImpExpr_sp_upper,	epi_ImpExpr_ppv_lower,	epi_ImpExpr_ppv_upper,	epi_ImpExpr_npv_lower,	epi_ImpExpr_npv_upper,	epi_ImpExpr_plr_lower,	epi_ImpExpr_plr_upper,	epi_ImpExpr_nlr_lower,	epi_ImpExpr_nlr_upper)
  
  list_epi_fisher_testresult = list()
  list_epi_fisher_testresult[["epi_fisher_testresult_colname"]] = epi_fisher_testresult_colname
  list_epi_fisher_testresult[["epi_fisher_testresult"]] = epi_fisher_testresult
  
  return(list_epi_fisher_testresult)
}



Get_com_PCsimpleUni<-function(PCsimple_mat,ImpExpr_mat)
{
  epi.tests_ImpExpr = epi.tests(ImpExpr_mat)
  epi_ImpExpr_aprev = epi.tests_ImpExpr$detail[1,2] # apparent prevalence.
  epi_ImpExpr_tprev = epi.tests_ImpExpr$detail[2,2]  # true prevalence.
  epi_ImpExpr_se = epi.tests_ImpExpr$detail[3,2]   # Sensitivity
  epi_ImpExpr_sp = epi.tests_ImpExpr$detail[4,2]   # Specificity
  epi_ImpExpr_ppv = epi.tests_ImpExpr$detail[9,2]   # Positive predictive value
  epi_ImpExpr_npv = epi.tests_ImpExpr$detail[10,2]   # Negative predictive value
  epi_ImpExpr_plr = epi.tests_ImpExpr$detail[11,2]   # Positive likelihood ratio
  epi_ImpExpr_nlr = epi.tests_ImpExpr$detail[12,2]   # Negative likelihood ratio
  
  ImpExpr_res = c(epi_ImpExpr_aprev,epi_ImpExpr_tprev,epi_ImpExpr_se,epi_ImpExpr_sp,epi_ImpExpr_ppv,epi_ImpExpr_npv,epi_ImpExpr_plr,epi_ImpExpr_nlr)
  ImpExpr_res_cols = c("epi_ImpExpr_aprev","epi_ImpExpr_tprev","epi_ImpExpr_se","epi_ImpExpr_sp",
                       "epi_ImpExpr_ppv","epi_ImpExpr_npv","epi_ImpExpr_plr","epi_ImpExpr_nlr")
  
  PCsimple_res = PCsimple_mat[["PCsimple_rval_combine"]]
  PCsimple_res_cols = PCsimple_mat[["PCsimple_rval_combine_colname"]]
  
  comb_res = c(ImpExpr_res,PCsimple_res)
  comb_res_cols = c(ImpExpr_res_cols,PCsimple_res_cols)
  
  list_epi_res_comp = list()
  list_epi_res_comp[["comb_res"]] = comb_res
  list_epi_res_comp[["comb_res_cols"]] = comb_res_cols
  return(list_epi_res_comp)
  
}



Get_related_variable_pvalue_opt<-function(dt_IDA,dt_IDA2,dt_IDA5)
{
  
  #save correlation above
  Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final,  use='complete.obs' ) 
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final)
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  #add info from pvalue optimizaiton for cov
  Vx2 = grep("true_jointIDA", colnames(dt_IDA2))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue  = cor(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue = RMSE(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue = cor(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue = RMSE(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  
  
  Vx2 = grep("true_jointIDA", colnames(dt_IDA5))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue_iterate  = cor(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue_iterate = RMSE(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate = cor(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate = RMSE(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  
  related_variable_pvalue_opt = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Pearson_pvalue,Cor_TureIDA_EstIDA_Pearson_pvalue_iterate,Cor_TureJointIDA_EstJointIDA_Pearson,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate ,RMSE_TrueIDA_EstIDA,RMSE_TrueIDA_EstIDA_pvalue,RMSE_TrueIDA_EstIDA_pvalue_iterate,RMSE_TrueJointIDA_EstJointIDA,RMSE_TrueJointIDA_EstJointIDA_pvalue,RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate)
  related_variable_pvalue_opt_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Pearson_pvalue","Cor_TureIDA_EstIDA_Pearson_pvalue_iterate","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_EstIDA_pvalue","RMSE_TrueIDA_EstIDA_pvalue_iterate","RMSE_TrueJointIDA_EstJointIDA","RMSE_TrueJointIDA_EstJointIDA_pvalue","RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate")
  
  list_related_variable_pvalue_opt = list()
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt"]] = related_variable_pvalue_opt
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt_colname"]] = related_variable_pvalue_opt_colname
  
  return(list_related_variable_pvalue_opt)
}


###Added by fyn. 2023.9.15. To be consistent with the changes of dt_IDA
Get_related_variable_pvalue_opt2<-function(dt_IDA,dt_IDA2,dt_IDA5)
{
  ### To format the data as 
  dt_IDA$true_IDA = as.numeric(dt_IDA$true_IDA)
  dt_IDA$Estimate_IDA = as.numeric(dt_IDA$Estimate_IDA)
  dt_IDA$true_jointIDA = as.numeric(dt_IDA$true_jointIDA)
  dt_IDA$est_jointIDA = as.numeric(dt_IDA$est_jointIDA)
  
  dt_IDA2$true_IDA = as.numeric(dt_IDA2$true_IDA)
  dt_IDA2$Estimate_IDA = as.numeric(dt_IDA2$Estimate_IDA)
  dt_IDA2$true_jointIDA = as.numeric(dt_IDA2$true_jointIDA)
  dt_IDA2$est_jointIDA = as.numeric(dt_IDA2$est_jointIDA)
  
  dt_IDA5$true_IDA = as.numeric(dt_IDA5$true_IDA)
  dt_IDA5$Estimate_IDA = as.numeric(dt_IDA5$Estimate_IDA)
  dt_IDA5$true_jointIDA = as.numeric(dt_IDA5$true_jointIDA)
  dt_IDA5$est_jointIDA = as.numeric(dt_IDA5$est_jointIDA)
  
  
  
  #save correlation above
  Cor_TureIDA_EstIDA_Pearson = cor(as.numeric(dt_IDA$true_IDA), 
                                   as.numeric(dt_IDA$Estimate_IDA),
                                   use='complete.obs') 
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(as.numeric(dt_IDA$true_jointIDA),
                                             as.numeric(dt_IDA$est_jointIDA), 
                                             use='complete.obs')
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA)
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  #add info from pvalue optimizaiton for cov
  Vx2 = grep("true_jointIDA", colnames(dt_IDA2))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue  = cor(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue = RMSE(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue = cor(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue = RMSE(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  
  
  Vx2 = grep("true_jointIDA", colnames(dt_IDA5))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue_iterate  = cor(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue_iterate = RMSE(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate = cor(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate = RMSE(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  
  related_variable_pvalue_opt = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Pearson_pvalue,Cor_TureIDA_EstIDA_Pearson_pvalue_iterate,Cor_TureJointIDA_EstJointIDA_Pearson,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate ,RMSE_TrueIDA_EstIDA,RMSE_TrueIDA_EstIDA_pvalue,RMSE_TrueIDA_EstIDA_pvalue_iterate,RMSE_TrueJointIDA_EstJointIDA,RMSE_TrueJointIDA_EstJointIDA_pvalue,RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate)
  related_variable_pvalue_opt_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Pearson_pvalue","Cor_TureIDA_EstIDA_Pearson_pvalue_iterate","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_EstIDA_pvalue","RMSE_TrueIDA_EstIDA_pvalue_iterate","RMSE_TrueJointIDA_EstJointIDA","RMSE_TrueJointIDA_EstJointIDA_pvalue","RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate")
  
  list_related_variable_pvalue_opt = list()
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt"]] = related_variable_pvalue_opt
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt_colname"]] = related_variable_pvalue_opt_colname
  
  return(list_related_variable_pvalue_opt)
}




#*****************************************************************************************
# Revised by yinly, Sep 26, 2023
# Description: evaluate the performance of multivariate linear regression 
#*****************************************************************************************
# @adj: true adjacency matrix
# @gene_index: gene index 
# @n3: phenotype index
# @ImpExpr_MLR_mat: coefficients of genes from multivariate linear regression with 5 cloumns: 1st(gene index), 2nd(coefficient),3rd(SE),4th(t value),5th(pvalue)
# @genes_num: number of genes
# @pval: pvalue cutoff for associated genes
#*****************************************************************************************
MVLR_Performance <- function(adj,gene_index,n3,ImpExpr_MLR_mat,genes_num,pval){
  Pheno_vec = adj[gene_index,n3]
  
  MVLR_sig_index <- which(ImpExpr_MLR_mat[,5]<pval) # The 5th column is the pvalue column
  MVLR_gene_index <- ImpExpr_MLR_mat[MVLR_sig_index,1] # extract detected susceptibility genes from MVLR 
  real_directCausal = gene_index[Pheno_vec!=0]
  
  cell_11 = length(intersect(MVLR_gene_index, real_directCausal))
  cell_12 = length(MVLR_gene_index) - cell_11
  cell_21 = length(real_directCausal) - cell_11
  cell_22 = genes_num - cell_11 - cell_12 - cell_21  
  
  MVLR_eval_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  MVLR_rval <- epi.tests(MVLR_eval_mat, conf.level = 0.95)
  MVLR_rval_aprev = MVLR_rval$detail[1,2] # apparent prevalence.
  MVLR_rval_tprev = MVLR_rval$detail[2,2]  # true prevalence.
  MVLR_rval_se = MVLR_rval$detail[3,2]   # Sensitivity
  MVLR_rval_sp = MVLR_rval$detail[4,2]   # Specificity
  MVLR_rval_ppv = MVLR_rval$detail[9,2]   # Positive predictive value
  MVLR_rval_npv = MVLR_rval$detail[10,2]   # Negative predictive value
  MVLR_rval_plr = MVLR_rval$detail[11,2]   # Positive likelihood ratio
  MVLR_rval_nlr = MVLR_rval$detail[12,2]   # Negative likelihood ratio
  
  MVLR_rval_combine = c(MVLR_rval_aprev,MVLR_rval_tprev,MVLR_rval_se,MVLR_rval_sp,MVLR_rval_ppv,MVLR_rval_npv,MVLR_rval_plr,MVLR_rval_nlr)
  MVLR_rval_combine_colname = c("MVLR_rval_aprev","MVLR_rval_tprev","MVLR_rval_se","MVLR_rval_sp","MVLR_rval_ppv","MVLR_rval_npv","MVLR_rval_plr","MVLR_rval_nlr")
  
  names(MVLR_rval_combine) <- MVLR_rval_combine_colname
  return(MVLR_rval_combine)   
}


#*************************************************************************************************
# Revised by yinly, Sep 26, 2023
# Description: evaluate the performance of multivariate linear regression with lasso/elastic-net
#*************************************************************************************************
# @adj: true adjacency matrix
# @gene_index: gene index 
# @n3: phenotype index
# @ImpExpr_Lasso_mat: coefficients of genes from mvlr_with_glmnet with 5 cloumns: 1st(gene index), 2nd(coefficient),3rd(SE),
#  4th(t value),5th(pvalue);coefficients of genes from mvlr_with_glmnet2 with 2 cloumns:1st(gene index), 2nd(coefficient)
# @genes_num: number of genes
# @target_column: define the coefficients column used for performance evaluation
#*************************************************************************************************
glmnet_Performance <- function(adj,gene_index,n3,ImpExpr_Lasso_mat,genes_num,target_column){
  Pheno_vec = adj[gene_index,n3]
  
  glmnet_sig_index <- which(ImpExpr_Lasso_mat[,target_column]!=0) # non-zero coefficient indicates association
  glmnet_gene_index <- ImpExpr_Lasso_mat[glmnet_sig_index,1] # extract detected susceptibility genes from glmnet 
  real_directCausal = gene_index[Pheno_vec!=0]
  
  cell_11 = length(intersect(glmnet_gene_index, real_directCausal))
  cell_12 = length(glmnet_gene_index) - cell_11
  cell_21 = length(real_directCausal) - cell_11
  cell_22 = genes_num - cell_11 - cell_12 - cell_21  
  
  glmnet_eval_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  glmnet_rval <- epi.tests(glmnet_eval_mat, conf.level = 0.95)
  glmnet_rval_aprev = glmnet_rval$detail[1,2] # apparent prevalence.
  glmnet_rval_tprev = glmnet_rval$detail[2,2]  # true prevalence.
  glmnet_rval_se = glmnet_rval$detail[3,2]   # Sensitivity
  glmnet_rval_sp = glmnet_rval$detail[4,2]   # Specificity
  glmnet_rval_ppv = glmnet_rval$detail[9,2]   # Positive predictive value
  glmnet_rval_npv = glmnet_rval$detail[10,2]   # Negative predictive value
  glmnet_rval_plr = glmnet_rval$detail[11,2]   # Positive likelihood ratio
  glmnet_rval_nlr = glmnet_rval$detail[12,2]   # Negative likelihood ratio
  
  glmnet_rval_combine = c(glmnet_rval_aprev,glmnet_rval_tprev,glmnet_rval_se,glmnet_rval_sp,glmnet_rval_ppv,glmnet_rval_npv,glmnet_rval_plr,glmnet_rval_nlr)
  glmnet_rval_combine_colname = c("glmnet_rval_aprev","glmnet_rval_tprev","glmnet_rval_se","glmnet_rval_sp","glmnet_rval_ppv","glmnet_rval_npv","glmnet_rval_plr","glmnet_rval_nlr")
  
  names(glmnet_rval_combine) <- glmnet_rval_combine_colname
  return(glmnet_rval_combine)   
}

