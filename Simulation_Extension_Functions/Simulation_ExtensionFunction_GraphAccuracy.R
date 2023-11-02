library(epiR)
#***************************************************************************************************
# Aim of this function is to perform epi.test of our overall causal graph 
#***************************************************************************************************
# @ - adj_true: adjacency matrix of true graph 
# @ - adj_est: adjacency matrix of estimated graph
# @ - beta4Fscore: settings of beta for caculating Fscore
#***************************************************************************************************
Cacul_Graph_Accuracy <- function(adj_true,adj_est,beta4Fscore = 0.5) {
  adj_minus = adj_true - adj_est
  adj2 = ifelse(adj_minus==0,0,8) # This code is used for distinguish the situation of C11 and c00. Added by fyn. 2023.8.25. 
  adj3 = adj_true + adj2
  
  # ref to the note in this link https://rdrr.io/cran/epiR/man/epi.tests.html
  c11 = sum(adj3==1) # There is an edge in both true and est graph.
  c01 = sum(adj_minus==-1) # There is an edge in est graph, but no edge in true graph.
  c10 = sum(adj_minus==1) # There is an edge in true graph, but no edge in est graph.
  c00 = sum(adj3==0)-ncol(adj_true) # There is no edge both in true and est graph. Added by fyn. 2023.10.4. But we need to exclude the nodes in the diagonal of the adj matrix
  
  # con_table = c(c11,c01,c10,c00)
  # con_table = table(con_table)
  dat <- as.table(matrix(c(c11,c01,c10,c00), nrow = 2, byrow = TRUE))
  colnames(dat) <- c("Dis+","Dis-")
  rownames(dat) <- c("Test+","Test-")
  
  
  # Need to check how to extract the value in graph_rval.  updated by fyn. 2023.8.25
  graph_rval <- epi.tests(dat, conf.level = 0.95)
  # Check the refs to see the format of the output of epi.test
  
  graph_rval_aprev = graph_rval$detail[1,2] # apparent prevalence.
  graph_rval_tprev = graph_rval$detail[2,2]  # true prevalence.
  graph_rval_se = graph_rval$detail[3,2]   # Sensitivity (recall)
  graph_rval_sp = graph_rval$detail[4,2]   # Specificity
  graph_rval_ppv = graph_rval$detail[9,2]   # Positive predictive value (precision)
  graph_rval_npv = graph_rval$detail[10,2]   # Negative predictive value
  graph_rval_plr = graph_rval$detail[11,2]   # Positive likelihood ratio
  graph_rval_nlr = graph_rval$detail[12,2]   # Negative likelihood ratio
  
  # Caculate the F1 score based on the precision and recall 
  F1_score = 2 * (graph_rval_ppv * graph_rval_se) / (graph_rval_ppv + graph_rval_se)
  tmp = beta4Fscore*beta4Fscore*graph_rval_ppv+graph_rval_se
  Fbeta_score = (1+beta4Fscore*beta4Fscore)*(graph_rval_ppv * graph_rval_se)/tmp
  
  F2_score = (1+2*2)*(graph_rval_ppv * graph_rval_se)/(2*2*graph_rval_ppv+graph_rval_se)
  F0.5_score = (1+0.5*0.5)*(graph_rval_ppv * graph_rval_se)/(0.5*0.5*graph_rval_ppv+graph_rval_se)
  
  
  graph_rval_combine = c(graph_rval_aprev,graph_rval_tprev,graph_rval_se,graph_rval_sp,graph_rval_ppv,graph_rval_npv,graph_rval_plr,graph_rval_nlr,
                         F1_score,F0.5_score,F2_score,Fbeta_score)
  graph_rval_combine_colname = c("graph_rval_aprev","graph_rval_tprev","graph_rval_se","graph_rval_sp","graph_rval_ppv","graph_rval_npv","graph_rval_plr","graph_rval_nlr",
                                 "F1_score","F0.5_score","F2_score","Fbeta_score")
  
  
  
  list_graph_rval_combine = list()
  list_graph_rval_combine[["graph_rval_combine"]] = graph_rval_combine
  list_graph_rval_combine[["graph_rval_combine_colname"]] = graph_rval_combine_colname
  return(list_graph_rval_combine)
}



#***************************************************************************************************
# This is the previous version of 'Cacul_Graph_Accuracy' function, with the edges in the diagonal of matrix
#***************************************************************************************************
# @ - adj_true: adjacency matrix of true graph 
# @ - adj_est: adjacency matrix of estimated graph
# @ - beta4Fscore: settings of beta for caculating Fscore
#***************************************************************************************************
Cacul_Graph_Accuracy_pre <- function(adj_true,adj_est,beta4Fscore = 0.5) {
  adj_minus = adj_true - adj_est
  adj2 = ifelse(adj_minus==0,0,8) # This code is used for distinguish the situation of C11 and c00. Added by fyn. 2023.8.25. 
  adj3 = adj_true + adj2
  
  # ref to the note in this link https://rdrr.io/cran/epiR/man/epi.tests.html
  c11 = sum(adj3==1) # There is an edge in both true and est graph.
  c01 = sum(adj_minus==-1) # There is an edge in est graph, but no edge in true graph.
  c10 = sum(adj_minus==1) # There is an edge in true graph, but no edge in est graph.
  c00 = sum(adj3==0) # There is no edge both in true and est graph.
  
  # con_table = c(c11,c01,c10,c00)
  # con_table = table(con_table)
  dat <- as.table(matrix(c(c11,c01,c10,c00), nrow = 2, byrow = TRUE))
  colnames(dat) <- c("Dis+","Dis-")
  rownames(dat) <- c("Test+","Test-")
  
  
  # Need to check how to extract the value in graph_rval.  updated by fyn. 2023.8.25
  graph_rval <- epi.tests(dat, conf.level = 0.95)
  # Check the refs to see the format of the output of epi.test
  
  graph_rval_aprev = graph_rval$detail[1,2] # apparent prevalence.
  graph_rval_tprev = graph_rval$detail[2,2]  # true prevalence.
  graph_rval_se = graph_rval$detail[3,2]   # Sensitivity (recall)
  graph_rval_sp = graph_rval$detail[4,2]   # Specificity
  graph_rval_ppv = graph_rval$detail[9,2]   # Positive predictive value (precision)
  graph_rval_npv = graph_rval$detail[10,2]   # Negative predictive value
  graph_rval_plr = graph_rval$detail[11,2]   # Positive likelihood ratio
  graph_rval_nlr = graph_rval$detail[12,2]   # Negative likelihood ratio
  
  # Caculate the F1 score based on the precision and recall 
  F1_score = 2 * (graph_rval_ppv * graph_rval_se) / (graph_rval_ppv + graph_rval_se)
  tmp = beta4Fscore*beta4Fscore*graph_rval_ppv+graph_rval_se
  Fbeta_score = (1+beta4Fscore*beta4Fscore)*(graph_rval_ppv * graph_rval_se)/tmp
  
  F2_score = (1+2*2)*(graph_rval_ppv * graph_rval_se)/(2*2*graph_rval_ppv+graph_rval_se)
  F0.5_score = (1+0.5*0.5)*(graph_rval_ppv * graph_rval_se)/(0.5*0.5*graph_rval_ppv+graph_rval_se)
  
  
  graph_rval_combine = c(graph_rval_aprev,graph_rval_tprev,graph_rval_se,graph_rval_sp,graph_rval_ppv,graph_rval_npv,graph_rval_plr,graph_rval_nlr,
                         F1_score,F0.5_score,F2_score,Fbeta_score)
  graph_rval_combine_colname = c("graph_rval_aprev","graph_rval_tprev","graph_rval_se","graph_rval_sp","graph_rval_ppv","graph_rval_npv","graph_rval_plr","graph_rval_nlr",
                                 "F1_score","F0.5_score","F2_score","Fbeta_score")
  
  
  
  list_graph_rval_combine = list()
  list_graph_rval_combine[["graph_rval_combine"]] = graph_rval_combine
  list_graph_rval_combine[["graph_rval_combine_colname"]] = graph_rval_combine_colname
  return(list_graph_rval_combine)
}



