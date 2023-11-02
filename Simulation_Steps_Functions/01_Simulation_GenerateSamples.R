#***************************************************************************************************
# Aim of this function is to simulate the gene-pheno causal network, which is marked as true graph in our paper
# The generated graph includes only one trait
#***************************************************************************************************
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
# @ - node_num: number of nodes in the simulated causal network
# @ - all_serial: sequence of 1 to number of nodes
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
#***************************************************************************************************
Simulate_Network<-function(s,g,p,node_num,all_serial,gene_index,snp_index) #prop_direct_causal_genes=prop_direct_causal_genes)
{
  
  rDAG <- randomDAG(n= node_num, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
  adj = as(rDAG, "matrix") #adj[i,j]!=0 represent that there is an edge i->j
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  rDAG = igraph.to.graphNEL(rDAG)
  
  tmp_loop = g* prop_direct_causal_genes 
  
  for (lp in 1:tmp_loop)
  {
    tmp_x = sample(gene_index, size = 1,replace = FALSE)
    
    tmp_y = n3 #only one phenotype
    
    adj_tmp = as(rDAG, "matrix") 
    
    weight_tmp = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[tmp_x,tmp_y] = weight_tmp
    
  }
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  
  rDAG = igraph.to.graphNEL(rDAG)
  
  non_zero_location = which(adj!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  
  #pre-processing in data
  #set SNP-SNP relation to 0
  #set gene-SNPs relation to 0
  
  for (i in 1:rows)
  {
    #genes and snps are randomly selected from data instead of choosing them in order as before
    
    
    #set SNP-SNP relation to 0
    if((non_zero_location[i,1] %in% snp_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    #set gene-SNPs relation to 0
    else if((non_zero_location[i,1] %in% gene_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    
  }
  
  #set pheno-other relation to 0
  adj[n3,]=0
  
  adj[snp_index,n3] = 0  #set snps to phenotype relation to 0
  
  
  #ensure that there is at least one snp related to gene.If there is a gene which it does not have related snps, the imputed gene expression data will be the same in all samples.
  for(gene in gene_index)
  {
    if(length(which(adj[snp_index,gene]!=0))==0)
    {
      snp_loc = sample(snp_index, 1, replace = FALSE)
      weight_tmp = runif(1, min = randomDAG_lB, max = randomDAG_uB)
      adj[snp_loc,gene] = weight_tmp
    }
  }
  
  
  #set snp-gene weaker to be consistent with real situation
  times_to_snpGeneWeight = randomDAG_uB/max_weight_snp_to_gene
  adj[snp_index,gene_index] = adj[snp_index,gene_index]/times_to_snpGeneWeight
  
  # adjust adj to amatType format, to make sure graph is DAG after update
  non_zero_location2 = which(adj!=0,arr.ind = T)
  
  #temp comment: edge with weight to be 1
  adj2 = adj
  for(i in 1:dim(non_zero_location2)[1])
  {
    a = non_zero_location2[i,1] # no-zero location in row
    b = non_zero_location2[i,2] # no-zero location in column
    adj2[a,b] = 1
  }
  
  amat = t(adj2) 
  # in amatType, a[i,j] = 0; a[j,i]=1, represent an edge i->j
  # reference https://www.rdocumentation.org/packages/pcalg/versions/2.6-2/topics/amatType
  
  isDAGresult = isValidGraph(amat = amat, type = "dag",verbose = TRUE)
  isCPDAGresult = isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG   completed partially directed acyclic graph 
  isPDAGresult = isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG
  
  if(isDAGresult||isCPDAGresult)
  {
    return(adj) #if the graph is DAG or CPDAG, then return.
  }
}



#***************************************************************************************************
# Aim of this function is to simulate the gene-pheno causal network, which is marked as true graph in our paper
# The generated graph includes both trait1 and trait2 and there is a directed causal relation between trait1 to trait2
#***************************************************************************************************
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
# @ - node_num: number of nodes in the simulated causal network
# @ - all_serial: sequence of 1 to number of nodes
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
#***************************************************************************************************
Simulate_Network_forMultiTraits<-function(s,g,p,node_num,all_serial,gene_index,snp_index) 
{
  
  rDAG <- randomDAG(n= node_num, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
  adj = as(rDAG, "matrix") #adj[i,j]!=0 represent that there is an edge i->j
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  rDAG = igraph.to.graphNEL(rDAG)
  
  tmp_loop = g* prop_direct_causal_genes 
  
  for (lp in 1:tmp_loop)
  {
    tmp_x = sample(gene_index, size = 1,replace = FALSE)
    
    tmp_y = n3 #only one phenotype
    
    adj_tmp = as(rDAG, "matrix") 
    
    weight_tmp = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[tmp_x,tmp_y] = weight_tmp
    
  }
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  
  rDAG = igraph.to.graphNEL(rDAG)
  
  non_zero_location = which(adj!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  
  #pre-processing in data
  #set SNP-SNP relation to 0
  #set gene-SNPs relation to 0
  
  for (i in 1:rows)
  {
    #genes and snps are randomly selected from data instead of choosing them in order as before
    #set SNP-SNP relation to 0
    if((non_zero_location[i,1] %in% snp_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    #set gene-SNPs relation to 0
    else if((non_zero_location[i,1] %in% gene_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    
  }
  
  #set pheno-other relation to 0
  adj[n3,]=0
  
  adj[snp_index,n3] = 0  #set snps to phenotype relation to 0
  
  
  #ensure that there is at least one snp related to gene.If there is a gene which it does not have related snps, the imputed gene expression data will be the same in all samples.
  for(gene in gene_index)
  {
    if(length(which(adj[snp_index,gene]!=0))==0)
    {
      snp_loc = sample(snp_index, 1, replace = FALSE)
      weight_tmp = runif(1, min = randomDAG_lB, max = randomDAG_uB)
      adj[snp_loc,gene] = weight_tmp
    }
  }
  
  
  #set snp-gene weaker to be consistent with real situation
  times_to_snpGeneWeight = randomDAG_uB/max_weight_snp_to_gene
  adj[snp_index,gene_index] = adj[snp_index,gene_index]/times_to_snpGeneWeight
  
  # adjust adj to amatType format, to make sure graph is DAG after update
  non_zero_location2 = which(adj!=0,arr.ind = T)
  
  ### when the number of phenotype is equal to 2
  index1 = node_num-1
  index2 = node_num
  ### Added by fyn 2023.9.24, in order to ensure that trait1 cause trait2
  adj[index1,index1] = 0
  adj[index1,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
  adj[index2,index1] = 0
  adj[index2,index2] = 0
  
  
  ### Added by fyn 2023.9.25, to test the problem of Genetic pleiotropy, ensure there is at least one gene which is associated with two trait.
  geneT1 = adj[gene_index,index1]
  geneT2 = adj[gene_index,index2]
  
  if(intersect(geneT1,geneT2)==0){
    random_i = sample(gene_index, 1, replace = FALSE)
    
    adj[random_i,index1] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[random_i,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    
  }
  
  
  #temp comment: edge with weight to be 1
  adj2 = adj
  adj2 = ifelse(adj2[]!=0,1,0) 
  
  
  # for(i in 1:dim(non_zero_location2)[1])
  # {
  #   a = non_zero_location2[i,1] # no-zero location in row
  #   b = non_zero_location2[i,2] # no-zero location in column
  #   adj2[a,b] = 1
  # }
  
  amat = t(adj2) 
  # in amatType, a[i,j] = 0; a[j,i]=1, represent an edge i->j
  # reference https://www.rdocumentation.org/packages/pcalg/versions/2.6-2/topics/amatType
  
  isDAGresult = isValidGraph(amat = amat, type = "dag",verbose = TRUE)
  isCPDAGresult = isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG   completed partially directed acyclic graph 
  isPDAGresult = isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG
  
  if(isDAGresult||isCPDAGresult)
  {
    return(adj) #if the graph is DAG or CPDAG, then return.
  }
  
}





#***************************************************************************************************
# Aim of this function is to simulate the gene-pheno causal network, which is marked as true graph in our paper
# Simulate multi traits and ensure a direct casual relation from trait1 to trait2, and there are some SNPs which can be selected as instrument variable
#***************************************************************************************************
# @ - adj: adjacency matrix of a random DAG
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
# @ - node_num: number of nodes in the simulated causal network
#***************************************************************************************************
Simulate_Network_forHP<-function(adj,gene_index,snp_index,s,g,p,node_num)
{
  rDAG <- randomDAG(n= node_num, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
  adj = as(rDAG, "matrix") #adj[i,j]!=0 represent that there is an edge i->j
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  rDAG = igraph.to.graphNEL(rDAG)
  
  tmp_loop = g* prop_direct_causal_genes 
  
  for (lp in 1:tmp_loop)
  {
    tmp_x = sample(gene_index, size = 1,replace = FALSE)
    
    tmp_y = n3 #only one phenotype
    
    adj_tmp = as(rDAG, "matrix") 
    
    weight_tmp = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[tmp_x,tmp_y] = weight_tmp
    
  }
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  
  rDAG = igraph.to.graphNEL(rDAG)
  
  non_zero_location = which(adj!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  
  #pre-processing in data
  #set SNP-SNP relation to 0
  #set gene-SNPs relation to 0
  
  for (i in 1:rows)
  {
    #genes and snps are randomly selected from data instead of choosing them in order as before
    #set SNP-SNP relation to 0
    if((non_zero_location[i,1] %in% snp_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    #set gene-SNPs relation to 0
    else if((non_zero_location[i,1] %in% gene_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    
  }
  
  #set pheno-other relation to 0
  adj[n3,]=0
  
  adj[snp_index,n3] = 0  #set snps to phenotype relation to 0
  
  
  #ensure that there is at least one snp related to gene.If there is a gene which it does not have related snps, the imputed gene expression data will be the same in all samples.
  for(gene in gene_index)
  {
    if(length(which(adj[snp_index,gene]!=0))==0)
    {
      snp_loc = sample(snp_index, 1, replace = FALSE)
      weight_tmp = runif(1, min = randomDAG_lB, max = randomDAG_uB)
      adj[snp_loc,gene] = weight_tmp
    }
  }
  
  
  #set snp-gene weaker to be consistent with real situation
  times_to_snpGeneWeight = randomDAG_uB/max_weight_snp_to_gene
  adj[snp_index,gene_index] = adj[snp_index,gene_index]/times_to_snpGeneWeight
  
  # adjust adj to amatType format, to make sure graph is DAG after update
  non_zero_location2 = which(adj!=0,arr.ind = T)
  
  ### when the number of phenotype is equal to 2
  index1 = node_num-1
  index2 = node_num
  ### Added by fyn 2023.9.24, in order to ensure that trait1 cause trait2
  adj[index1,index1] = 0
  adj[index1,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
  adj[index2,index1] = 0
  adj[index2,index2] = 0
  
  
  ### Added by fyn 2023.9.25, to test the problem of Genetic pleiotropy, ensure there is at least one gene which is associated with two trait.
  geneT1 = adj[gene_index,index1]
  geneT2 = adj[gene_index,index2]
  
  if(intersect(geneT1,geneT2)==0){
    random_i = sample(gene_index, 1, replace = FALSE)
    
    adj[random_i,index1] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[random_i,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
  }
  
  
  ### Added by fyn 2023.9.30. To ensure we have SNP which can be selected as instrument variable when apply MR method
  SNPT1 = adj[snp_index,index1]
  random_snp_index = sample(snp_index, 1, replace = FALSE)
  adj[random_snp_index,index1] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
  
  geneT2_index = adj[gene_index,index2]
  nonZero_index = geneT2_index[geneT2_index!=0]
  gene_sel_index = sample(nonZero_index, 1, replace = FALSE)
  gene_sel_index2 = as.numeric(names(gene_sel_index))
  ### Added the edge from SNP to geneX, and geneX cause T2
  adj[random_snp_index,gene_sel_index2] = runif(1, min = randomDAG_lB, max = randomDAG_uB)
  
  
  #temp comment: edge with weight to be 1
  adj2 = adj
  adj2 = ifelse(adj2[]!=0,1,0) 
  
  
  # for(i in 1:dim(non_zero_location2)[1])
  # {
  #   a = non_zero_location2[i,1] # no-zero location in row
  #   b = non_zero_location2[i,2] # no-zero location in column
  #   adj2[a,b] = 1
  # }
  
  amat = t(adj2) 
  # in amatType, a[i,j] = 0; a[j,i]=1, represent an edge i->j
  # reference https://www.rdocumentation.org/packages/pcalg/versions/2.6-2/topics/amatType
  
  isDAGresult = isValidGraph(amat = amat, type = "dag",verbose = TRUE)
  isCPDAGresult = isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG   completed partially directed acyclic graph 
  isPDAGresult = isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG
  
  if(isDAGresult||isCPDAGresult)
  {
    adj_info_list = list()
    adj_info_list[["adj"]] = adj
    adj_info_list[["random_snp_index"]] = random_snp_index
    adj_info_list[["gene_sel_index2"]] = gene_sel_index2
    return(adj_info_list) #if the graph is DAG or CPDAG, then return.
  }
  
}





#***************************************************************************************************
# Aim of this function is to simulate the gene-pheno causal network, which is marked as true graph in our paper
# Increase the weight of SNP to outcome and genes, this function is more used for testing, not in the main simulated framework
#***************************************************************************************************
# @ - adj: adjacency matrix of a random DAG
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of phenotypes in the simulated causal network
# @ - node_num: number of nodes in the simulated causal network
#***************************************************************************************************
Simulate_Network_forMultiTraits_SNPwIncrease<-function(s,g,p,node_num,all_serial,gene_index,snp_index) 
{
  
  rDAG <- randomDAG(n= node_num, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
  adj = as(rDAG, "matrix") #adj[i,j]!=0 represent that there is an edge i->j
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  rDAG = igraph.to.graphNEL(rDAG)
  
  tmp_loop = g* prop_direct_causal_genes 
  
  for (lp in 1:tmp_loop)
  {
    tmp_x = sample(gene_index, size = 1,replace = FALSE)
    
    tmp_y = n3 #only one phenotype
    
    adj_tmp = as(rDAG, "matrix") 
    
    weight_tmp = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[tmp_x,tmp_y] = weight_tmp
    
  }
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  
  rDAG = igraph.to.graphNEL(rDAG)
  
  non_zero_location = which(adj!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  
  #pre-processing in data
  #set SNP-SNP relation to 0
  #set gene-SNPs relation to 0
  
  for (i in 1:rows)
  {
    #genes and snps are randomly selected from data instead of choosing them in order as before
    #set SNP-SNP relation to 0
    if((non_zero_location[i,1] %in% snp_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    #set gene-SNPs relation to 0
    else if((non_zero_location[i,1] %in% gene_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    
  }
  
  #set pheno-other relation to 0
  adj[n3,]=0
  
  adj[snp_index,n3] = 0  #set snps to phenotype relation to 0
  
  
  #ensure that there is at least one snp related to gene.If there is a gene which it does not have related snps, the imputed gene expression data will be the same in all samples.
  for(gene in gene_index)
  {
    if(length(which(adj[snp_index,gene]!=0))==0)
    {
      snp_loc = sample(snp_index, 1, replace = FALSE)
      weight_tmp = runif(1, min = randomDAG_lB, max = randomDAG_uB)
      adj[snp_loc,gene] = weight_tmp
    }
  }
  
  ### Updated by fyn. 2023.10.4. We do NOT set the snp-gene to be weaker
  # #set snp-gene weaker to be consistent with real situation
  # times_to_snpGeneWeight = randomDAG_uB/max_weight_snp_to_gene
  # adj[snp_index,gene_index] = adj[snp_index,gene_index]/times_to_snpGeneWeight
  
  # adjust adj to amatType format, to make sure graph is DAG after update
  non_zero_location2 = which(adj!=0,arr.ind = T)
  
  ### when the number of phenotype is equal to 2
  index1 = node_num-1
  index2 = node_num
  ### Added by fyn 2023.9.24, in order to ensure that trait1 cause trait2
  adj[index1,index1] = 0
  adj[index1,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
  adj[index2,index1] = 0
  adj[index2,index2] = 0
  
  
  ### Added by fyn 2023.9.25, to test the problem of Genetic pleiotropy, ensure there is at least one gene which is associated with two trait.
  geneT1 = adj[gene_index,index1]
  geneT2 = adj[gene_index,index2]
  
  if(intersect(geneT1,geneT2)==0){
    random_i = sample(gene_index, 1, replace = FALSE)
    
    adj[random_i,index1] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[random_i,index2] = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    
  }
  
  
  #temp comment: edge with weight to be 1
  adj2 = adj
  adj2 = ifelse(adj2[]!=0,1,0) 
  
  
  # for(i in 1:dim(non_zero_location2)[1])
  # {
  #   a = non_zero_location2[i,1] # no-zero location in row
  #   b = non_zero_location2[i,2] # no-zero location in column
  #   adj2[a,b] = 1
  # }
  
  amat = t(adj2) 
  # in amatType, a[i,j] = 0; a[j,i]=1, represent an edge i->j
  # reference https://www.rdocumentation.org/packages/pcalg/versions/2.6-2/topics/amatType
  
  isDAGresult = isValidGraph(amat = amat, type = "dag",verbose = TRUE)
  isCPDAGresult = isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG   completed partially directed acyclic graph 
  isPDAGresult = isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG
  
  if(isDAGresult||isCPDAGresult)
  {
    return(adj) #if the graph is DAG or CPDAG, then return.
  }
  
}



#***************************************************************************************************
# Aim of this function is to simulate noise when simulate the causal networks without hidden variable
#***************************************************************************************************
# @ - sample_num2: number of samples when generating the causal networks
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
#***************************************************************************************************
Simulate_Noise<-function(sample_num2,s,g)
{
  Sigma <- diag(s) #snps error : construct a unit diagonal matrix
  eMat_s <- mvrnorm(sample_num2, mu = rep(0, s), Sigma = Sigma) #produce snp error
  
  eMat_g = NULL
  tmp_num = g
  
  for (g_num in 1:tmp_num)
  {
    eMat_g_tmp = rnorm(sample_num2, mean = 0, sd = eMat_g_tmp_rnorm_sd) # set variance of gene data error manually
    eMat_g =cbind(eMat_g,eMat_g_tmp)
  }
  
  eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
  eMat = cbind(eMat_s,eMat_g,eMat_p)
  
  eMat[,gene_index] = eMat_g
  eMat[,node_num] = eMat_p
  eMat[,snp_index]=eMat_s
  
  return(eMat)
}


#***************************************************************************************************
# Aim of this function is to simulate noise when simulate the causal networks with hidden variable
#***************************************************************************************************
# @ - sample_num2: number of samples when generating the causal networks
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - node_num: number of nodes in the simulated causal network
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
#***************************************************************************************************
Simulate_Noise_HiddenVar<-function(sample_num2,s,g,node_num,gene_index,snp_index)
{
  Sigma <- diag(s) #snps error : construct a unit diagonal matrix
  eMat_s <- mvrnorm(sample_num2, mu = rep(0, s), Sigma = Sigma) #produce snp error
  
  eMat_g = NULL
  tmp_num = g
  
  for (g_num in 1:tmp_num)
  {
    eMat_g_tmp = rnorm(sample_num2, mean = 0, sd = eMat_g_tmp_rnorm_sd) # set variance of gene data error manually
    eMat_g =cbind(eMat_g,eMat_g_tmp)
  }
  
  eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
  eMat = cbind(eMat_s,eMat_g,eMat_p)
  
  eMat[,gene_index] = eMat_g
  eMat[,node_num] = eMat_p
  eMat[,snp_index]=eMat_s
  
  return(eMat)
}



#***************************************************************************************************
# Aim of this function is to simulate noise when simulate the causal networks with multiple traits in the graph
#***************************************************************************************************
# @ - sample_num2: number of samples when generating the causal networks
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p_num: number of traits(or phenotypes) in the simulated causal network
# @ - gene_index: index of genes in the simulated causal network
# @ - snp_index : index of snps in the simulated causal network
#***************************************************************************************************
Simulate_Noise_ForMultiTraits<-function(sample_num2,s,g,p_num,snp_index,gene_index)
{
  Sigma <- diag(s) #snps error : construct a unit diagonal matrix
  eMat_s <- mvrnorm(sample_num2, mu = rep(0, s), Sigma = Sigma) #produce snp error
  
  eMat_g = NULL
  tmp_num = g
  
  for (g_num in 1:tmp_num)
  {
    eMat_g_tmp = rnorm(sample_num2, mean = 0, sd = eMat_g_tmp_rnorm_sd) # set variance of gene data error manually
    eMat_g =cbind(eMat_g,eMat_g_tmp)
  }
  
  eMat_all_p = NULL
  for(p_index in 1:p_num)
  {
    eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
    eMat_all_p = cbind(eMat_p,eMat_all_p)
  }
  # colnames(eMat_all_p) = paste0("eMat_p",seq(1,p_num))
  
  
  # eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
  eMat = cbind(eMat_s,eMat_g,eMat_all_p)
  total_num = s+g+p_num
  colnames(eMat) = seq(1,total_num)
  
  
  eMat[,snp_index]=eMat_s
  eMat[,gene_index] = eMat_g
  # eMat[,node_num] = eMat_p
  
  return(eMat)
}




#***************************************************************************************************
# Aim of this function is to simulate noise when simulate the causal networks 
# With settings of different parameter of standard deviation of gene expression data
#***************************************************************************************************
# @ - sample_num2: number of samples when generating the causal networks
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - eMat_g_tmp_rnorm_sd: standard deviation of gene expression data
#***************************************************************************************************
Simulate_Noise_DiffError<-function(sample_num2,s,g,eMat_g_tmp_rnorm_sd)
{
  Sigma <- diag(s) #snps error : construct a unit diagonal matrix
  eMat_s <- mvrnorm(sample_num2, mu = rep(0, s), Sigma = Sigma) #produce snp error
  
  eMat_g = NULL
  tmp_num = g
  
  for (g_num in 1:tmp_num)
  {
    eMat_g_tmp = rnorm(sample_num2, mean = 0, sd = eMat_g_tmp_rnorm_sd) # set variance of gene data error manually
    eMat_g =cbind(eMat_g,eMat_g_tmp)
  }
  
  eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
  eMat = cbind(eMat_s,eMat_g,eMat_p)
  
  eMat[,gene_index] = eMat_g
  eMat[,node_num] = eMat_p
  eMat[,snp_index]=eMat_s
  
  return(eMat)
}



#***************************************************************************************************
# Aim of this function is to simulate noise when simulate the causal networks 
# simulate different measurement error in case and control groups
#***************************************************************************************************
# @ - sample_num2: number of samples when generating the causal networks
# @ - s: number of SNPs in the simulated causal network
# @ - g: number of genes in the simulated causal network
# @ - p: number of traits(or phenotypes) in the simulated causal network
# @ - gene_index: gene_index in the simulated causal network 
# @ - gene_expression_error: the error needed to simulate the differential and non differential measurement error
#***************************************************************************************************
Adjust_Error_CaseControl<-function(samples_data2,gene_expression_error,s,g,p,gene_index)
{
  outcome_index = ncol(samples_data2)
  colnames(samples_data2)[ncol(samples_data2)] =  "outcome"# The last column should be the value of the outcome
  
  outcome_value = samples_data2[,ncol(samples_data2)] 
  mean_outcome = mean(outcome_value)
  samples_data2 = as.data.table(samples_data2)
  samples_data2$casecontrol_tag = ifelse(samples_data2$outcome>=mean_outcome,1,0) # Using mean as the cutoff to distinguish case and control
  
  
  
  samples_data2 = as.matrix(samples_data2)
  # tmp = samples_data2[,genei]
  # > head(tmp)
  # [1] -2.8974951  0.4004909 -2.9154572 -0.7136560 -1.6077103 -0.4807393
  
  ### Added by fyn. 2023.10.2. Add more error to case group
  for(genei in gene_index)
  {
    samples_data2[,genei] = ifelse(samples_data2[,"casecontrol_tag"]==1,
                                   (samples_data2[,genei]+gene_expression_error),
                                   samples_data2[,genei])
    # > head(samples_data2[,genei])
    # [1] -2.8974951  0.4004909 -2.9154572 -0.6636560 -1.6077103 -0.4807393
    # > head(samples_data2[,"casecontrol_tag"])
    # [1] 0 0 0 1 0 0
  }
  
  return(samples_data2)
  
}

