library(data.table)

args <- commandArgs(trailingOnly = TRUE)
directory = args[1]
run = args[2]
is.SE = args[3]

#########  read in input files #################
consolidated_GLM_file = paste(directory,run,"_GLM_outputs_consolidated.txt",sep="")
consolidated_GLM = fread(consolidated_GLM_file,sep="\t",header=TRUE)

if (is.SE==1){
  SICILIAN_junctions_file = paste(directory,run,".tsv",sep="")
  SICILIAN_junctions = fread(SICILIAN_junctions_file,sep="\t",header = TRUE)
  SICILIAN_junctions[,c("frac_passed_samples","emp.p_glmnet_constrained_median","called"):=NULL]
}
###############################################

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
  consolidated_GLM = merge(SICILIAN_junctions,consolidated_GLM_uniq[,list(refName_newR1,frac_passed_samples,emp.p_glmnet_constrained_median,called)],by.x="refName_newR1",by.y="refName_newR1",all.x=TRUE,all.y=FALSE)
} else{
  setnames(consolidated_GLM,"sample_name","cell")
  consolidated_GLM[,channel:=cell]
  consolidated_GLM[,called:=0]
  consolidated_GLM[(emp.p_glmnet_corrected_constrained_median<0.1) & ((frac_passed_samples>0.5) | (frac_passed_samples==0.5 & sample_num>2)),called:=1]
  
  consolidated_GLM[, chrR1A:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][1], by = refName_newR1]
  consolidated_GLM[, chrR1B:=strsplit(refName_newR1, split = "[:|]")[[1]][5], by = refName_newR1]
  consolidated_GLM[, juncPosR1A:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][3]), by = refName_newR1]
  consolidated_GLM[, juncPosR1B:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][6]), by = refName_newR1]
  consolidated_GLM[, geneR1A_uniq:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][2], by = refName_newR1]
  
  options(scipen = 999)
  round.bin = 100000
  consolidated_GLM[geneR1A_uniq=="",geneR1A_uniq:="unknown"]
  consolidated_GLM[geneR1A_uniq=="unknown", binR1A:=round(juncPosR1A/round.bin)*round.bin]
  consolidated_GLM[geneR1A_uniq=="unknown",geneR1A_uniq:=paste("unknown",chrR1A,binR1A,sep = "_")]
  consolidated_GLM[, binR1A:= NULL]
  consolidated_GLM[, binR1B:= NULL]
}

write.table(consolidated_GLM,paste(directory,"/",run,"_with_postprocessing.txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
