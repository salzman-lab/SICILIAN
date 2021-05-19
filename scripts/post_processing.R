library(data.table)

args <- commandArgs(trailingOnly = TRUE)
directory = args[1]
run = args[2] 

consolidated_GLM_file = paste(directory,run,"_GLM_outputs_consolidated.txt",sep="")
SICILIAN_junctions_file = paste(directory,run,".tsv",sep="")
is.SE = 1

#########  read in input files #################
consolidated_GLM = fread(consolidated_GLM_file,sep="\t",header=TRUE)
SICILIAN_junctions = fread(SICILIAN_junctions_file,sep="\t",header = TRUE)
###############################################
print(names(consolidated_GLM))
SICILIAN_junctions[,c("frac_passed_samples","emp.p_glmnet_constrained_median","called"):=NULL]

#SICILIAN_junctions = data.table()

#list_files = c("Krasnow_COVID_10x_pilot4_201214S2.tsv","Krasnow_COVID_10x_pilot4_201214S3.tsv","Krasnow_COVID_10x_pilot4_201214S4.tsv","Krasnow_COVID_10x_pilot4_201214S5.tsv","Krasnow_COVID_10x_pilot4_201214S6.tsv","Krasnow_COVID_10x_pilot4_201214S7.tsv","Krasnow_COVID_10x_pilot4_201214S8.tsv","Krasnow_COVID_10x_pilot4_201214S9.tsv","Krasnow_COVID_10x_pilot4_201214S10.tsv")

#for (counter in 1: length(list_files)){
# file_name = list_files[counter]
# report = fread(paste(SICILIAN_junctions_file,file_name, sep = ""), sep = "\t",header=TRUE)
# SICILIAN_junctions = rbind(SICILIAN_junctions,report)
#}

consolidated_GLM[,passed_tot:=0]
consolidated_GLM[(postprocess_passed==1) & (frac_genomic_reads<0.1),passed_tot:=1]
consolidated_GLM[,frac_passed_samples:=mean(passed_tot),by=refName_newR1]
consolidated_GLM[,frac_passed_samples:=max(frac_passed_samples),by=refName_newR1]
consolidated_GLM[,max_numread:=max(numReads),by=refName_newR1]
consolidated_GLM[max_numread==1,frac_passed_samples:=0]
consolidated_GLM[,sample_num:=length(unique(sample_name)),by=refName_newR1]

# consolidated_GLM_uniq has the median emp.p and the fraction of passed samples based on postprocessing for each unique junction in the list
if (is.SE==1){
  consolidated_GLM_uniq = consolidated_GLM[!duplicated(refName_newR1),list(refName_newR1,frac_passed_samples,sample_num,emp.p_glmnet_constrained_median)]
  consolidated_GLM_uniq[,called:=0]
  consolidated_GLM_uniq[(emp.p_glmnet_constrained_median<0.1) & ((frac_passed_samples>0.5) | (frac_passed_samples==0.5 & sample_num>2)),called:=1]
  consolidated_GLM_uniq = consolidated_GLM_uniq[!duplicated(refName_newR1)]
  SICILIAN_junctions = merge(SICILIAN_junctions,consolidated_GLM_uniq[,list(refName_newR1,frac_passed_samples,emp.p_glmnet_constrained_median,called)],by.x="refName_newR1",by.y="refName_newR1",all.x=TRUE,all.y=FALSE)
} else{
  consolidated_GLM_uniq = consolidated_GLM[!duplicated(refName_newR1),list(refName_newR1,frac_passed_samples,sample_num,emp.p_glmnet_corrected_constrained_median)]
  consolidated_GLM_uniq[,called:=0]
  consolidated_GLM_uniq[(emp.p_glmnet_corrected_constrained_median<0.1) & ((frac_passed_samples>0.5) | (frac_passed_samples==0.5 & sample_num>2)),called:=1]
  consolidated_GLM_uniq = consolidated_GLM_uniq[!duplicated(refName_newR1)]
  SICILIAN_junctions = merge(SICILIAN_junctions,consolidated_GLM_uniq[,list(refName_newR1,frac_passed_samples,emp.p_glmnet_corrected_constrained_median,called)],by.x="refName_newR1",by.y="refName_newR1",all.x=TRUE,all.y=FALSE)
}

write.table(SICILIAN_junctions,paste(directory,"/",run,"_with_postprocessing.txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
