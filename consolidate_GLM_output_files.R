#!/usr/bin/env Rscript
#This script consolidates all the final report files for knife or machete across a TCGA cancer or a GTEx tissue into a single file that will be later used in the machete_SBT_post_process.R script to incorporate SBT detection frequency
## consolidate_machete_results_all_samples.R "spork" "GTEx" "Lung"
library(data.table)
library(stringr)
library(cutpointr)
library(ggplot2)
library(glmnet)
library(dplyr)

postprocess_model_SS <- function(consolidated_GLM){
  ## selecting training data
  consolidated_GLM[( (frac_multimapping>0.8) | (frac_genomic_reads>0.8) | (frac_anomaly==1)) & (numReads>1),train:=0]
  consolidated_GLM[(frac_multimapping==0) &(frac_genomic_reads==0) & (frac_anomaly==0)  & (!refName_newR1%like%"chrM") &(numReads>1),train:=1]
  
  n.neg = nrow(consolidated_GLM[train==0])
  n.pos = nrow(consolidated_GLM[train==1])
  n.neg = min(n.neg,200000)
  n.pos = min(n.pos,200000)  # number of positive reads that we want to subsample from the list of all reads
  all_neg_reads = which(consolidated_GLM$train==0)
  all_pos_reads = which(consolidated_GLM$train==1)
  consolidated_GLM[,train:=NULL]
  consolidated_GLM[sample(all_neg_reads, n.neg, replace= FALSE), train := 0]
  consolidated_GLM[sample(all_pos_reads, n.pos, replace= FALSE), train := 1]
  
  class.weight = min(n.pos, n.neg)
  consolidated_GLM[train == 0, cur_weight := n.pos/n.neg]
  consolidated_GLM[train == 1, cur_weight := 1]
  
  
  # compute the refined junction noise score
  round.bin = 50
  consolidated_GLM[, juncPosR1A:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][3]), by = refName_newR1]
  consolidated_GLM[, juncPosR1B:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][6]), by = refName_newR1]
  consolidated_GLM[, chrR1A:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][1], by = refName_newR1]
  consolidated_GLM[, chrR1B:=strsplit(refName_newR1, split = "[:|]")[[1]][5], by = refName_newR1]
  consolidated_GLM[, binR1A:=round(juncPosR1A/round.bin)*round.bin]
  consolidated_GLM[, binR1B:=round(juncPosR1B/round.bin)*round.bin]
  consolidated_GLM[,bin_chrR1A:=paste(binR1A, chrR1A),by=refName_newR1]
  consolidated_GLM[,bin_chrR1B:=paste(binR1B, chrR1B),by=refName_newR1]
  consolidated_GLM[fileTypeR1=="Aligned", njunc_binR1A_new:=length(unique(refName_newR1)), by = bin_chrR1A]
  consolidated_GLM[fileTypeR1=="Aligned", njunc_binR1B_new:=length(unique(refName_newR1)), by = bin_chrR1B]
  consolidated_GLM[,c("bin_chrR1A","bin_chrR1B"):=NULL]
  consolidated_GLM[,max_njunc:=max(njunc_binR1A_new,njunc_binR1B_new),by=1:nrow(consolidated_GLM)]
  
  
  ## now we train and apply the postprocessing model
  consolidated_GLM[,discrete_overlap:=round((median_overlap_R1)/5)]
  regression_formula = as.formula("train ~  discrete_overlap*ave_entropyR1+ave_max_run_R1   + max_njunc")
  x_glmnet = model.matrix(regression_formula, consolidated_GLM[!is.na(train)])
  glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(consolidated_GLM[!is.na(train)]$train), family =c("binomial"),weights = consolidated_GLM[!is.na(train)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits=c(Inf,Inf,Inf,0,0,Inf), lower.limits=c(-Inf,0,0,-Inf,-Inf,0) )
  coefficients = coef(glmnet_model_constrained)
  print(coefficients)
  f_glmnet = model.matrix(update(regression_formula, refName_newR1 ~ .), consolidated_GLM)
  a=predict(object = glmnet_model_constrained,newx = f_glmnet,type = "response", s = "lambda.1se", se.fit = TRUE)
  consolidated_GLM$postprocess_prob=a
  cp = cutpointr(consolidated_GLM[!is.na(train)], postprocess_prob, train,method = maximize_metric, metric = sum_sens_spec)
  consolidated_GLM[postprocess_prob<cp$optimal_cutpoint,postprocess_passed:=0]
  consolidated_GLM[postprocess_prob>=cp$optimal_cutpoint,postprocess_passed:=1]
  consolidated_GLM[ave_max_run_R1 > 14,postprocess_passed:=0]
  consolidated_GLM[ave_entropyR1 < 3,postprocess_passed:=0]
  if(length(which(!coefficients==0)) < 2){
    consolidated_GLM[,postprocess_passed:=1]
    consolidated_GLM[median_overlap_R1<10,postprocess_passed:=0]
    consolidated_GLM[ave_max_run_R1 > 9,postprocess_passed:=0]
    consolidated_GLM[ave_entropyR1 < 4,postprocess_passed:=0]
  }
  
  if (is.SE==0){
    consolidated_GLM[frac_anomaly == 1,postprocess_passed:=0]
  }
  consolidated_GLM[frac_multimapping == 1,postprocess_passed:=0]
  
  consolidated_GLM[,c("binR1A","binR1B","njunc_binR1A_new","njunc_binR1B_new","cur_weight","train"):=NULL]
  return(consolidated_GLM)
}


args <- commandArgs(trailingOnly = TRUE)
directory = args[1]
run = args[2]     #TS_pilot_10X_withinbam          
exon_pickle = args[3]
splice_pickle = args[4]
is.SE = args[5]   #run = "TS_pilot_10X_withinbam"
#################### Inputs (usually should not be changed unles the format of the file locations have changed##########################
#directory = paste("/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/",run,"/",sep = "")
consolidated_list_name = paste(run,"GLM_outputs_consolidated.txt",sep = "_")
###################################################

list_files = list.files(directory, pattern = "GLM_output.txt", all.files = TRUE,recursive = TRUE)  # the list of all GLM output files
consolidated_list = data.table()


for (counter in 1: length(list_files)){
  
  file_name = list_files[[counter]]
  sample_name =  strsplit(file_name,split = "/",fixed = TRUE)[[1]][1]   #the name of the sample will be added as a column
  print(file_name)
  if(is.SE==1){
    #"frac_multimapping"
    report = fread(paste(directory,file_name, sep = ""), sep = "\t", header = TRUE, fill = TRUE, select = c("refName_newR1","postprocess_passed","frac_genomic_reads","numReads","median_overlap_R1","sd_overlap","ave_max_run_R1","ave_entropyR1","emp.p_glmnet_constrained","fileTypeR1"))
    report = report[fileTypeR1=="Aligned"]
  } else{
    report = fread(paste(directory,file_name, sep = ""), sep = "\t", header = TRUE, fill = TRUE, select = c("refName_newR1","frac_genomic_reads","frac_multimapping","frac_anomaly","numReads","median_overlap_R1","sd_overlap","ave_max_run_R1","ave_entropyR1","emp.p_glmnet_corrected_constrained","fileTypeR1"))
    report = report[fileTypeR1=="Aligned"]
  }
  report[,sample_name := sample_name]
  #  report = report[,list(refName_newR1, sample_name, numReads, p_predicted_glm, junc_cdf_glm, p_predicted_glmnet, junc_cdf_glmnet, p_predicted_glmnet_constrained, junc_cdf_glmnet_constrained, emp.p_glm, emp.p_glmnet, emp.p_glmnet_constrained, seqR1)]
  consolidated_list = rbind(consolidated_list,report)
}

if (is.SE==1){
  consolidated_list[,emp.p_glmnet_constrained_median:=median(emp.p_glmnet_constrained),by=refName_newR1]
} else {
  consolidated_list[,emp.p_glmnet_corrected_constrained_median:= median(emp.p_glmnet_corrected_constrained),by=refName_newR1]
}


if (is.SE==1){
  write.table(consolidated_list,paste(directory,consolidated_list_name,sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
} else{
  consolidated_list = postprocess_model_SS(consolidated_list) # here I perform the postprocessing model for the SS2 across the entire consolidated GLM file
  write.table(consolidated_list,paste(directory,consolidated_list_name,sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
  script_directory = getwd()
  system(paste("python3 ",script_directory,"/scripts/ann_splices.py -i ",directory,consolidated_list_name, " -o ",directory,consolidated_list_name," -e ",exon_pickle, " -s ", splice_pickle,sep = "")) # here I annotate junctions for the SS2 data
}

