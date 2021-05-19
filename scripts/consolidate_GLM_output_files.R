#!/usr/bin/env Rscript
#This script consolidates all the final report files for knife or machete across a TCGA cancer or a GTEx tissue into a single file that will be later used in the machete_SBT_post_process.R script to incorporate SBT detection frequency
## consolidate_machete_results_all_samples.R "spork" "GTEx" "Lung"

require(data.table)
require(stringr)

args <- commandArgs(trailingOnly = TRUE)
directory = args[1]
run = args[2]     #TS_pilot_10X_withinbam          
is.SE = args[3]   #run = "TS_pilot_10X_withinbam"
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
  consolidated_list[,emp.p_glmnet_corrected_constrained_median:= median(emp.p_glmnet_constrained),by=refName_newR1]
}
write.table(consolidated_list,paste(directory,consolidated_list_name,sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
