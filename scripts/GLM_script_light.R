#!/usr/bin/env Rscript
## TheGLM script that first predict a per-read probability for each read alignment and then compute an aggregated score for each junction 

# load and install required packages for the script

if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("glmnet")) {
  install.packages("glmnet", dependencies = TRUE)
  library(glmnet)
}
if (!require("tictoc")) {
  install.packages("tictoc", dependencies = TRUE)
  library(tictoc)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

if (!require("cutpointr")) {
  install.packages("cutpointr", dependencies = TRUE)
  library(stringr)
}

if (!require("GenomicAlignments")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GenomicAlignments")
  
  library(GenomicAlignments)
}



tic("Entire code")

postprocessing_model <- function(GLM_output, class_input, is.SE, is.10X){
  # the postprocessing_model assigns a probability to each junction based on the junction-level features and that probability can be used for calling junctions
  
  GLM_output_chimeric = GLM_output[fileTypeR1=="Chimeric"] # the model is only for aligned junctions but will keep the chimeric junctions as well to add them at the end
  GLM_output = GLM_output[fileTypeR1=="Aligned"]
  
  class_input[NHR1A>1,multimapping:=1]
  class_input[NHR1A==1,multimapping:=0]
  class_input[,frac_multimapping:=mean(multimapping),by=refName_newR1]
  class_input_uniq=class_input[!duplicated(refName_newR1)]
  GLM_output[,frac_multimapping:=NULL]
  GLM_output[,frac_mutimapping:=NULL]
  GLM_output = merge(GLM_output,class_input_uniq[,list(refName_newR1,frac_multimapping)],all.x=TRUE,all.y=FALSE,by.x="refName_newR1",by.y="refName_newR1")
  GLM_output[,train:=NULL]
  if (is.SE==1){
    GLM_output[(frac_multimapping==0) &(frac_genomic_reads==0)  & (!refName_newR1%like%"chrM") &(numReads>1),train:=1]
    GLM_output[( (frac_multimapping>0.8) | (frac_genomic_reads>0.8)) & (numReads>1),train:=0]
  } else {
    GLM_output[(frac_multimapping==0) &(frac_genomic_reads==0) & (frac_anomaly==0)  & (!refName_newR1%like%"chrM") &(numReads>1),train:=1]
    GLM_output[( (frac_multimapping>0.8) | (frac_genomic_reads>0.8) | (frac_anomaly==1)) & (numReads>1),train:=0]
  }
  n.neg = nrow(GLM_output[train == 0])
  n.pos = nrow(GLM_output[train == 1])
  class.weight = min(n.pos, n.neg)
  GLM_output[train == 0, cur_weight := n.pos/n.neg]
  GLM_output[train == 1, cur_weight := 1]
  
  
  # compute the refined junction noise score
  round.bin = 50
  GLM_output[, binR1A:=round(juncPosR1A/round.bin)*round.bin]
  GLM_output[, binR1B:=round(juncPosR1B/round.bin)*round.bin]
  GLM_output[,bin_chrR1A:=paste(binR1A, chrR1A)]
  GLM_output[,bin_chrR1B:=paste(binR1B, chrR1B)]
  GLM_output[fileTypeR1=="Aligned", njunc_binR1A_new:=length(unique(refName_newR1)), by = bin_chrR1A]
  GLM_output[fileTypeR1=="Aligned", njunc_binR1B_new:=length(unique(refName_newR1)), by = bin_chrR1B]
  GLM_output[,c("bin_chrR1A","bin_chrR1B"):=NULL]
  GLM_output[,max_njunc:=max(njunc_binR1A_new,njunc_binR1B_new),by=1:nrow(GLM_output)]
  
  if(nrow(GLM_output[train==0])<30 | nrow(GLM_output[train==1])<30){
    GLM_output[,postprocess_passed:=1]
    GLM_output[median_overlap_R1<10,postprocess_passed:=0]
    GLM_output[ave_max_run_R1 > 9,postprocess_passed:=0]
    GLM_output[ave_entropyR1 < 4,postprocess_passed:=0]
    if (is.SE==0){
      GLM_output[frac_anomaly == 1,postprocess_passed:=0]
    }
    GLM_output[frac_multimapping == 1,postprocess_passed:=0]
  } else{
    ## now we train and apply the postprocessing model
    GLM_output[,discrete_overlap:=round((median_overlap_R1)/5)]
    regression_formula = as.formula("train ~  discrete_overlap*ave_entropyR1+ave_max_run_R1   + max_njunc")
    x_glmnet = model.matrix(regression_formula, GLM_output[!is.na(train)])
    glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(GLM_output[!is.na(train)]$train), family =c("binomial"),weights = GLM_output[!is.na(train)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits=c(Inf,Inf,Inf,0,0,Inf), lower.limits=c(-Inf,0,0,-Inf,-Inf,0) )
    coefficients = coef(glmnet_model_constrained)
    print(coefficients) 
    f_glmnet = model.matrix(update(regression_formula, refName_newR1 ~ .), GLM_output)
    a=predict(object = glmnet_model_constrained,newx = f_glmnet,type = "response", s = "lambda.1se", se.fit = TRUE)
    GLM_output$postprocess_prob=a
    cp = cutpointr(GLM_output[!is.na(train)], postprocess_prob, train,method = maximize_metric, metric = sum_sens_spec)
    GLM_output[postprocess_prob<cp$optimal_cutpoint,postprocess_passed:=0]
    GLM_output[postprocess_prob>=cp$optimal_cutpoint,postprocess_passed:=1]
    coefficients = as.vector(coefficients)
    if(length(which(coefficients!=0)) < 2){
      GLM_output[,postprocess_passed:=1]
      GLM_output[median_overlap_R1<10,postprocess_passed:=0]
      GLM_output[ave_max_run_R1 > 9,postprocess_passed:=0]
      GLM_output[ave_entropyR1 < 4,postprocess_passed:=0]
    }
    if (is.SE==0){
      GLM_output[frac_anomaly == 1,postprocess_passed:=0]
    }
    GLM_output[frac_multimapping == 1,postprocess_passed:=0]
  }
  GLM_output[ave_max_run_R1 > 14,postprocess_passed:=0]
  GLM_output[ave_entropyR1 < 3,postprocess_passed:=0]
  
  GLM_output = dplyr::bind_rows(GLM_output_chimeric,GLM_output)
  return(GLM_output)
  
}


postprocessing_hardthreshold <- function(GLM_output, class_input, is.10X, is.SE){
  
  called_junctions = GLM_output[fileTypeR1=="Aligned"]
  called_junctions[,intron_length:=abs(juncPosR1A-juncPosR1B), by = 1:nrow(called_junctions)]
  
  #### filtering out the GLM report file for calling fusions
  if (is.SE == 0){
    called_junctions = called_junctions[(fileTypeR1=="Aligned") & numReads > 1 & emp.p_glmnet_corrected_constrained < 0.1 & frac_genomic_reads < 0.1]
  } else{
    called_junctions = called_junctions[(fileTypeR1=="Aligned") & numReads > 1 & emp.p_glmnet_constrained < 0.1 & frac_genomic_reads < 0.1]
  }
  called_junctions = called_junctions[sd_overlap > 0]
  called_junctions = called_junctions[intron_length > 30]
  called_junctions = called_junctions[ave_AT_run_R1 < 11]
  called_junctions = called_junctions[ave_entropyR1 > 3]
  
  if (is.10X == 1){
    class_input[,numReads_per_cell:=.N,by=paste(refName_newR1,barcode,sep="")]
    class_input = class_input[!duplicated(paste(refName_newR1,barcode,sep=""))]
    called_junctions = merge(class_input[,list(refName_newR1,barcode,numReads_per_cell)],called_junctions,by.x="refName_newR1",by.y="refName_newR1",all.x=FALSE,all.y=TRUE)
  }
  
  called_junctions = data.frame(called_junctions)
  
  if (is.SE==0){   # I need to have different column names here becuse of setcolorder
    cols_to_keep = c("refName_newR1","numReads","geneR1A_uniq","geneR1B_uniq","chrR1A","juncPosR1A","juncPosR1B","gene_strandR1A","intron_length","seqR1","sd_overlap","frac_genomic_reads","ave_AT_run_R1","ave_GC_run_R1","ave_max_run_R1","junc_cdf_glmnet_corrected_constrained","emp.p_glmnet_corrected_constrained","ave_entropyR1") 
  } else {
    cols_to_keep = c("refName_newR1","numReads","numReads_per_cell","barcode","geneR1A_uniq","geneR1B_uniq","chrR1A","juncPosR1A","juncPosR1B","gene_strandR1A","intron_length","seqR1","sd_overlap","frac_genomic_reads","ave_AT_run_R1","ave_GC_run_R1","ave_max_run_R1","junc_cdf_glmnet_constrained","emp.p_glmnet_constrained","ave_entropyR1") 
  }
  
  called_junctions = called_junctions[,which(names(called_junctions)%in%cols_to_keep)]
  called_junctions = setcolorder(called_junctions,cols_to_keep[cols_to_keep%in%names(called_junctions)])
  
  return(called_junctions)
}



#this function annotates missing/inserted domains for the called junctions
domain_annotation <- function(ucsc_domain_file,called_junctions){
  
  ucsc = fread(ucsc_domain_file, sep = "\t", header = FALSE)
  
  ###### extracting exonic positions in a domain  #######
  ucsc = ucsc[,list(V2,V3,V4,V5,V12,V13)]
  ucsc[,V3 := V3 + 1] #as I want that the start positions shows the first base in exon
  ucsc[,num_blocks := length(strsplit(V13,split = ",")[[1]]),by = V13]
  ucsc_for_introns = ucsc[num_blocks>1]  # I want this to obtain the intronic regions in the protein domain
  ucsc = data.frame(ucsc)
  ucsc = ucsc[rep(row.names(ucsc), ucsc$num_blocks), 1:7]
  ucsc$id = rownames(ucsc)
  ucsc = data.table(ucsc)
  ucsc[,exon_id := as.numeric(strsplit(id,split=".",fixed=TRUE)[[1]][2]),by=id]
  ucsc[,domain_id := as.numeric(strsplit(id,split=".",fixed=TRUE)[[1]][1]),by=id]
  ucsc[is.na(exon_id),exon_id := 0]
  ucsc[,exon_id := exon_id+1]
  ucsc[,domain_exon_start := V3 + as.numeric(strsplit(V13,split = ",",fixed = TRUE)[[1]][exon_id]),by = 1:nrow(ucsc)]
  ucsc[,domain_exon_end := V3 + as.numeric(strsplit(V13,split = ",",fixed = TRUE)[[1]][exon_id]) + as.numeric(strsplit(V12,split = ",",fixed = TRUE)[[1]][exon_id])-1,by = 1:nrow(ucsc)]
  ucsc_for_introns[,num_introns:=num_blocks-1]
  ucsc_for_introns = data.frame(ucsc_for_introns)
  ucsc_for_introns = ucsc_for_introns[rep(row.names(ucsc_for_introns), ucsc_for_introns$num_introns), 1:8]
  ucsc_for_introns$id = rownames(ucsc_for_introns)
  ucsc_for_introns = data.table(ucsc_for_introns)
  ucsc_for_introns[,intron_id := as.numeric(strsplit(id,split=".",fixed=TRUE)[[1]][2]),by=id]
  ucsc_for_introns[,domain_id := as.numeric(strsplit(id,split=".",fixed=TRUE)[[1]][1]),by=id]
  ucsc_for_introns[is.na(intron_id),intron_id := 0]
  ucsc_for_introns[,intron_id := intron_id+1]
  ucsc_for_introns[,domain_intron_start := V3 + as.numeric(strsplit(V13,split = ",",fixed = TRUE)[[1]][intron_id]) + as.numeric(strsplit(V12,split = ",",fixed = TRUE)[[1]][intron_id]),by = 1:nrow(ucsc_for_introns)]
  ucsc_for_introns[,domain_intron_end := V3 + as.numeric(strsplit(V13,split = ",",fixed = TRUE)[[1]][intron_id + 1]) -1, by = 1:nrow(ucsc_for_introns)]
  ucsc_for_introns = ucsc_for_introns[domain_intron_start<domain_intron_end]
  ###############################
  
  
  ucsc_extracted= ucsc[,list(V2,V5,domain_exon_start,domain_exon_end,domain_id,exon_id)]
  setnames(ucsc_extracted,old=c("V2","domain_exon_start","domain_exon_end"),new=c("chrR1A","juncPosR1A_new", "juncPosR1B_new"))
  setkey(ucsc_extracted, chrR1A, juncPosR1A_new, juncPosR1B_new)
  
  called_junctions = data.table(called_junctions)
  called_junctions_uniq = called_junctions[!duplicated(refName_newR1),list(chrR1A,juncPosR1A,juncPosR1B,refName_newR1)]
  called_junctions_uniq[,juncPosR1A_new:=min(juncPosR1A,juncPosR1B),by=1:nrow(called_junctions_uniq)]
  called_junctions_uniq[,juncPosR1B_new:=max(juncPosR1A,juncPosR1B),by=1:nrow(called_junctions_uniq)]
  called_junctions_uniq[,juncPosR1A_new:=juncPosR1A_new + 1] #Ineed the intronic coordinate as I want to overlap with domains
  called_junctions_uniq[,juncPosR1B_new:=juncPosR1B_new - 1] #Ineed the intronic coordinate as I want to overlap with domains
  
  
  out = foverlaps(called_junctions_uniq, ucsc_extracted, type="any")
  out[is.na(V5),V5:=""]
  out=out[!duplicated(paste(domain_id,refName_newR1))]
  out[,missing_domains:=paste(V5,collapse = ":"),by=refName_newR1]
  out = out[,list(refName_newR1, missing_domains)]
  out = out[!duplicated(refName_newR1)]
  called_junctions = merge(called_junctions,out,by.x="refName_newR1",by.y="refName_newR1",all.x=TRUE,all.y=TRUE)
  
  ## now we do the same overlaping process to find insertions in the domain
  called_junctions_uniq[,juncPosR1A_new:=min(juncPosR1A,juncPosR1B),by=1:nrow(called_junctions_uniq)]
  called_junctions_uniq[,juncPosR1B_new:=max(juncPosR1A,juncPosR1B),by=1:nrow(called_junctions_uniq)]
  called_junctions_uniq[,juncPosR1A_new_1:=juncPosR1A_new - 2]
  called_junctions_uniq[,juncPosR1B_new_1:=juncPosR1B_new + 2]
  ucsc_for_introns_extracted= ucsc_for_introns[,list(V2,V5,domain_intron_start,domain_intron_end,domain_id,intron_id)]
  
  # first we check the overlap for R1A side of the junction
  setnames(ucsc_for_introns_extracted,old=c("V2","domain_intron_start","domain_intron_end"),new=c("chrR1A","juncPosR1A_new_1", "juncPosR1A_new"))
  setkey(ucsc_for_introns_extracted, chrR1A, juncPosR1A_new_1, juncPosR1A_new)
  out_R1A = foverlaps(called_junctions_uniq, ucsc_for_introns_extracted[,list(chrR1A,juncPosR1A_new_1,juncPosR1A_new,V5,domain_id)], type="any")
  out_R1A = out_R1A[,list(refName_newR1,V5,domain_id)]
  
  # now we check the overlap for the R1B side of the junction
  setnames(ucsc_for_introns_extracted,old=c("chrR1A","juncPosR1A_new_1", "juncPosR1A_new"),new=c("chrR1A","juncPosR1B_new", "juncPosR1B_new_1"))
  setkey(ucsc_for_introns_extracted, chrR1A, juncPosR1B_new, juncPosR1B_new_1)
  out_R1B = foverlaps(called_junctions_uniq, ucsc_for_introns_extracted[,list(chrR1A,juncPosR1B_new,juncPosR1B_new_1,V5,domain_id)], type="any")
  out_R1B = out_R1B[,list(refName_newR1,V5,domain_id)]
  
  out_introns = rbind(out_R1A,out_R1B)
  out_introns[is.na(V5),V5:=""]
  out_introns=out_introns[!duplicated(paste(domain_id,refName_newR1))]
  out_introns[,num_inser_domains:=length(unique(domain_id)),by=refName_newR1]
  out_introns = out_introns[!((num_inser_domains>1) & (V5==""))]
  out_introns[,domain_insertions:=paste(V5,collapse = ":"),by=refName_newR1]
  out_introns = out_introns[,list(refName_newR1, domain_insertions)]
  out_introns = out_introns[!duplicated(refName_newR1)]
  
  called_junctions = merge(called_junctions,out_introns,by.x="refName_newR1",by.y="refName_newR1",all.x=TRUE,all.y=TRUE)
  return(called_junctions)
}


# this function is used for adding gene ensembl anmes and htseq counts from STAR output
add_ensembl <- function(gtf_file,directory,class_input,is.SE){
  
  gtf_info = fread(gtf_file,header = FALSE,sep="\t")
  gtf_entry_types=unique(gtf_info$V3)
  if ("gene" %in% gtf_entry_types){
    gtf_info[V3 == "gene",length := V5-V4]
  }
  
  gtf_info = gtf_info[(V3 == "exon") | (V3 == "gene")]
  
  # I want to find the mapping between gene id and gene name fir the genes in the gtf file
  gtf_info[,gene_id:=strsplit(strsplit(V9,split="gene_id ")[[1]][2],split=";")[[1]][1],by=V9]
  gtf_info[,gene_id:=gsub(" ","",gene_id), by = gene_id]
  gtf_info[,gene_id:=gsub("\"","",gene_id), by = gene_id]
  gtf_info = gtf_info[!duplicated(gene_id)]
  gtf_info[,gene_id := strsplit(gene_id,split = "." , fixed = TRUE)[[1]][1], by = gene_id]
  gtf_info[,gene_name := strsplit(strsplit(V9,split = "gene_name ")[[1]][2],split = ";")[[1]][1], by = V9]
  gtf_info[,gene_name := gsub("\"","",gene_name),by = gene_name]
  gene_name_id = gtf_info[,list(gene_name,gene_id)]
  
  if ("length" %in% names(gtf_info)){
    gene_length = gtf_info[V3 == "gene",list(gene_id,length)]
  }
  
  if(is.SE == 1){
    genecount_file = paste(directory,list.files(directory, pattern = "2ReadsPerGene.out.tab", all.files = FALSE),sep = "")
  } else {
    genecount_file = paste(directory,list.files(directory, pattern = "1ReadsPerGene.out.tab", all.files = FALSE),sep = "")
  }
  
  ######### read in gene count file #######
  gene_count = fread(genecount_file,sep = "\t",header = FALSE, skip = 4)
  if (gene_count[1,V1]%like%"ensembl"){
    gene_count = fread(genecount_file,sep = "\t",header = TRUE)
  }
  
  if (! ("ensembl_id" %in% names(gene_count)) ){  
    gene_count[,ensembl_id:=strsplit(V1,split = ".",fixed = TRUE)[[1]][1],by = 1:nrow(gene_count)]
  }
  if ("gene_name" %in% names(gene_count)){
    gene_count[,c("gene_name","length") := NULL]  # if there is a gene name column from the past in the file, we want to delete it to avoid errors
  }
  ##########################################
  
  
  
  #add HUGO gene names to the gene count file
  gene_count = merge(gene_count,gene_name_id,by.x = "ensembl_id",by.y = "gene_id",all.x = TRUE,all.y = FALSE)
  #gene_count[ensembl_id%like%"ERCC",gene_name:=ensembl_id]
  #gene_count[(ensembl_id %in% c("L","S","E","M","N")) | (ensembl_id %like% "ORF"),gene_name:=ensembl_id]
  gene_count[,V1 := NULL]
  gene_count = gene_count[!duplicated(gene_name)]
  gene_count = gene_count[!duplicated(ensembl_id)]
  total_read_unstranded = sum(gene_count$V2)
  total_read_stranded = sum(gene_count$V3)
  
  if("length" %in% names(gtf_info)){
    gene_count = merge(gene_count,gene_length,by.x = "ensembl_id",by.y = "gene_id",all.x = TRUE,all.y = FALSE)
    gene_count[,RPKM_unstranded:= V2/ (total_read_unstranded/(10^6) * length/(10^3))]
    gene_count[,RPKM_stranded:= V3/ (total_read_stranded/(10^6) * length/(10^3))]
  }
  
  # now add gene ensembl and gene counts to the class input file
  if ( "geneR1B_ensembl" %in% names(class_input) ){
    class_input[,geneR1A_name := NULL]
    class_input[,geneR1B_name := NULL]
    class_input[,geneR1A_ensembl := NULL]
    class_input[,geneR1B_ensembl := NULL]
    class_input[,geneR1A_expression_stranded := NULL] 
    class_input[,geneR1B_expression_stranded := NULL]
    class_input[,geneR1A_expression_unstranded := NULL]
    class_input[,geneR1B_expression_unstranded := NULL]
    class_input[,geneR1A_expression := NULL]
    class_input[,geneR1B_expression := NULL]
    class_input[,geneR1A_RPKM_stranded := NULL]
    class_input[,geneR1B_RPKM_stranded := NULL]
    class_input[,geneR1A_RPKM_unstranded := NULL]
    class_input[,geneR1B_RPKM_unstranded := NULL]
  }
  
  class_input = merge(class_input,unique(gene_name_id[!duplicated(gene_name),list(gene_name,gene_id)]),by.x = "geneR1A_uniq",by.y = "gene_name",all.x = TRUE,all.y = FALSE)
  setnames(class_input,old = "gene_id" ,new = "geneR1A_ensembl")
  class_input = merge(class_input,unique(gene_name_id[!duplicated(gene_name),list(gene_name,gene_id)]),by.x = "geneR1B_uniq",by.y = "gene_name",all.x = TRUE,all.y = FALSE)
  setnames(class_input,old = "gene_id" ,new = "geneR1B_ensembl")
  if("length" %in% names(gtf_info)){
    class_input = merge(class_input,gene_count[!duplicated(ensembl_id),list(ensembl_id,V2,V3,RPKM_unstranded,RPKM_stranded)],by.x = "geneR1A_ensembl",by.y = "ensembl_id",all.x = TRUE,all.y = FALSE)
    setnames(class_input,old = c("V2","V3","RPKM_unstranded","RPKM_stranded") ,new = c("geneR1A_expression_unstranded","geneR1A_expression_stranded","geneR1A_RPKM_unstranded","geneR1A_RPKM_stranded"))
    class_input = merge(class_input,gene_count[!duplicated(ensembl_id),list(ensembl_id,V2,V3,RPKM_unstranded,RPKM_stranded)],by.x = "geneR1B_ensembl",by.y = "ensembl_id",all.x = TRUE,all.y = FALSE)
    setnames(class_input,old = c("V2","V3","RPKM_unstranded","RPKM_stranded") ,new = c("geneR1B_expression_unstranded","geneR1B_expression_stranded","geneR1B_RPKM_unstranded","geneR1B_RPKM_stranded"))
  }else{
    class_input = merge(class_input,gene_count[!duplicated(ensembl_id),list(ensembl_id,V2,V3)],by.x = "geneR1A_ensembl",by.y = "ensembl_id",all.x = TRUE,all.y = FALSE)
    setnames(class_input,old = c("V2","V3") ,new = c("geneR1A_expression_unstranded","geneR1A_expression_stranded"))
    class_input = merge(class_input,gene_count[!duplicated(ensembl_id),list(ensembl_id,V2,V3)],by.x = "geneR1B_ensembl",by.y = "ensembl_id",all.x = TRUE,all.y = FALSE)
    setnames(class_input,old = c("V2","V3") ,new = c("geneR1B_expression_unstranded","geneR1B_expression_stranded"))
  }
  
  ## write output files
  write.table(gene_count,genecount_file,row.names = FALSE,sep = "\t",quote = FALSE)
  return(class_input)
}


compare_classinput_STARChimOut <- function(directory,is.SE){
  
  options(scipen = 999)  # this will make sure that the modified coordinates in chimeric or SJ files won't be written in scientific representation as we want to compare them with those in the class input file 
  star_SJ_output_1_file = list.files(directory,pattern = "2SJ.out.tab")
  star_chimeric_output_1_file = list.files(directory,pattern = "2Chimeric.out.junction")
  #  star_fusion_file = list.files(directory,pattern = "star-fusion.fusion_predictions.abridged.tsv" , recursive = TRUE)
  
  if (is.SE == 0){
    star_SJ_output_1_file = list.files(directory,pattern = "1SJ.out.tab")
    star_SJ_output_2_file = list.files(directory,pattern = "2SJ.out.tab")
    star_chimeric_output_1_file = list.files(directory,pattern = "1Chimeric.out.junction")
    star_chimeric_output_2_file = list.files(directory,pattern = "2Chimeric.out.junction")
  }
  
  ####### read in files ############
  chimeric1 =  fread(paste(directory,star_chimeric_output_1_file,sep = ""),sep = "\t",header = FALSE)
  star_SJ_output_1 =  fread(paste(directory,star_SJ_output_1_file,sep = ""),sep = "\t",header = FALSE)
  #  star_fusion = fread(paste(directory,star_fusion_file,sep = ""),sep = "\t" , header = TRUE)
  #############################
  
  # if the script has been run previously on the class inout file, delete the following columns to avoid duplicate column names
  class_input[,intron_motif:=NULL]
  class_input[,is.annotated:=NULL]
  class_input[,num_uniq_map_reads:=NULL]
  class_input[,num_multi_map_reads:=NULL]
  class_input[,maximum_SJ_overhang:=NULL]
  
  #converting intronic positions to exonic positions for chimeric junction output file for R1 
  chimeric1[V3 == "+", V2_exonic := V2 - 1]
  chimeric1[V3 == "-", V2_exonic := V2 + 1]
  chimeric1[V6 == "+", V5_exonic := V5 + 1]
  chimeric1[V6 == "-", V5_exonic := V5 - 1]
  chimeric1[,junction := paste(V1,V2_exonic,V4,V5_exonic,sep = ":")]
  
  star_SJ_output_1[,junction := paste(V1,V2,V1,V3,sep = ":")]
  
  class_input[,min_junc_pos:=min(juncPosR1A,juncPosR1B),by=paste(juncPosR1A,juncPosR1B)] # I do this to consistently have the minimum junc position first in the junction id used for comparing between the class input file and the STAR output file
  class_input[,max_junc_pos:=max(juncPosR1A,juncPosR1B),by=paste(juncPosR1A,juncPosR1B)]
  class_input[fileTypeR1 == "Chimeric",junction_compatible := paste(chrR1A,juncPosR1A,chrR1B,juncPosR1B,sep = ":")]
  class_input[fileTypeR1 == "Aligned",junction_compatible := paste(chrR1A,as.numeric(min_junc_pos) + 1,chrR1B,as.numeric(max_junc_pos) - 1,sep = ":")]
  class_input[,min_junc_pos:=NULL]
  class_input[,max_junc_pos:=NULL]
  
  class_input[fileTypeR1 == "Chimeric",is.STAR_Chim := 1]
  class_input[fileTypeR1 == "Chimeric" & !(junction_compatible %in% chimeric1$junction) ,is.STAR_Chim := 0]
  
  class_input[fileTypeR1 == "Aligned",is.STAR_SJ := 1]
  class_input[fileTypeR1 == "Aligned" & !(junction_compatible %in% star_SJ_output_1$junction) ,is.STAR_SJ := 0] 
  class_input = merge(class_input,star_SJ_output_1[,list(junction,V5,V6,V7,V8,V9)],by.x = "junction_compatible",by.y = "junction",all.x = TRUE,all.y = FALSE )
  
  class_input[fileTypeR1 == "Chimeric", is.STAR_Fusion := 0]
  
  in_star_chim_not_in_classinput = chimeric1[!(junction %in% class_input$junction_compatible)]
  in_star_SJ_not_in_classinput = star_SJ_output_1[ !(junction %in% class_input$junction_compatible)]
  
  class_input[,junction_compatible := NULL]
  setnames(class_input,old = c("V5","V6","V7","V8","V9"), new = c("intron_motif","is.annotated","num_uniq_map_reads","num_multi_map_reads","maximum_SJ_overhang"))
  ################################################
  
  #### write output files
  write.table(in_star_chim_not_in_classinput,paste(directory,"in_star_chim_not_in_classinput.txt",sep = ""),quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(in_star_SJ_not_in_classinput,paste(directory,"in_star_SJ_not_in_classinput.txt",sep = ""),quote = FALSE, row.names = FALSE, sep = "\t")
  return(class_input)
}



compute_class_error <- function(train_class, glm_predicted_prob){
  totalerr = sum(abs(train_class - round(glm_predicted_prob)))
  print (paste("total reads:", length(train_class)))
  print(paste("both negative", sum(abs(train_class+round(glm_predicted_prob))==0), "out of ", length(which(train_class==0))))
  print(paste("both positive", sum(abs(train_class+round(glm_predicted_prob))==2), "out of ", length(which(train_class==1))))
  print(paste("classification errors for glm", totalerr, "out of", length(train_class), totalerr/length(train_class) ))
}

compute_junc_cdf <- function(class_input, p_predicted_column, per_read_column, junc_cdf_column){
  # compute the junc_cdf scores
  
  
  names(class_input)[names(class_input)==per_read_column]="per_read_prob"
  names(class_input)[names(class_input)==p_predicted_column]="p_predicted"
  
  class_input[p_predicted==1, p_predicted :=0.999999999999]
  class_input[p_predicted==0, p_predicted :=10^-30]
  class_input[per_read_prob==1, per_read_prob:=0.999999999999]
  class_input[per_read_prob==0, per_read_prob:=10^-30]
  class_input[, log_per_read_prob:=log((1-per_read_prob) / per_read_prob), by = per_read_prob]
  setkey(class_input,refName_newR1)
  iter=10000
  class_input[, sum_log_per_read_prob:= sum(log_per_read_prob), by = refName_newR1]
  
  mu_i = mean(log( (1-class_input[fileTypeR1 == "Aligned"| (fileTypeR1 == "Chimeric" & (chrR1A == chrR1B) & (gene_strandR1A==gene_strandR1B) &  ((gene_strandR1A== "+" & juncPosR1A < juncPosR1B) | (gene_strandR1A== "-" & juncPosR1A > juncPosR1B)) & abs(juncPosR1A- juncPosR1B)<1000000 )]$per_read_prob)/ class_input[fileTypeR1 == "Aligned"| (fileTypeR1 == "Chimeric" & (chrR1A == chrR1B) & (gene_strandR1A==gene_strandR1B) &  ((gene_strandR1A== "+" & juncPosR1A < juncPosR1B) | (gene_strandR1A== "-" & juncPosR1A > juncPosR1B)) & abs(juncPosR1A- juncPosR1B)<1000000 )]$per_read_prob) )
  var_i = var(log( (1-class_input[fileTypeR1 == "Aligned"| (fileTypeR1 == "Chimeric" & (chrR1A == chrR1B) & (gene_strandR1A==gene_strandR1B) &  ((gene_strandR1A== "+" & juncPosR1A < juncPosR1B) | (gene_strandR1A== "-" & juncPosR1A > juncPosR1B)) & abs(juncPosR1A- juncPosR1B)<1000000 )]$per_read_prob)/ class_input[fileTypeR1 == "Aligned"| (fileTypeR1 == "Chimeric" & (chrR1A == chrR1B) & (gene_strandR1A==gene_strandR1B) &  ((gene_strandR1A== "+" & juncPosR1A < juncPosR1B) | (gene_strandR1A== "-" & juncPosR1A > juncPosR1B)) & abs(juncPosR1A- juncPosR1B)<1000000 )]$per_read_prob) )
  all_per_read_probs = class_input[(fileTypeR1 == "Aligned") | (fileTypeR1 == "Chimeric" & (chrR1A == chrR1B) & (gene_strandR1A==gene_strandR1B) &  ((gene_strandR1A== "+" & juncPosR1A < juncPosR1B) | (gene_strandR1A== "-" & juncPosR1A > juncPosR1B)) & abs(juncPosR1A- juncPosR1B)<1000000 )]$per_read_prob
  
  num_per_read_probs = length(all_per_read_probs)
  for (num_reads in 1:15){
    rnd_per_read_probs = matrix(0, iter, num_reads)
    rnd_per_read_probs = apply(rnd_per_read_probs,1, function(x) all_per_read_probs[sample(num_per_read_probs, num_reads)])
    rnd_per_read_probs = t(rnd_per_read_probs)
    if(num_reads == 1){
      rnd_per_read_probs = t(rnd_per_read_probs)  # for num_reads =1 I need to transepose twice since first I have a vector
    }
    null_dist = apply(rnd_per_read_probs,1, function(x) 1/( exp(sum(log( (1 - x)/x ))) + 1))
    null_dist[which(null_dist==1)]=0.999999999999
    null_dist[which(null_dist==0)]=10^-30
    class_input[fileTypeR1 == "Aligned" & numReads == num_reads, junc_cdf :=length(which(null_dist <= p_predicted))/iter, by = p_predicted]
  }
  
  class_input[fileTypeR1 == "Aligned" & numReads > 15, junc_cdf :=pnorm(sum_log_per_read_prob, mean = numReads*mu_i, sd = sqrt(numReads*var_i), lower.tail = FALSE), by = refName_newR1]
  
  names(class_input)[names(class_input) == "per_read_prob"] = per_read_column
  names(class_input)[names(class_input) == "p_predicted"] =  p_predicted_column
  names(class_input)[names(class_input) == "junc_cdf"] =  junc_cdf_column
  return(class_input)
}

uniformity_test <- function(dt, min_R1_offset, max_R1_offset) {
  possible_values = data.frame(vals = min_R1_offset:max_R1_offset)
  t=data.frame(table(dt))
  possible_values = merge(possible_values, t, by.x = "vals", by.y = "dt", all.x = TRUE, all.y = FALSE)
  possible_values = data.table(possible_values)
  possible_values[is.na(Freq)]$Freq = 0
  test = chisq.test(possible_values$Freq, simulate.p.value = TRUE,B = 200)
  #  test = chisq.test(possible_values$Freq)
  return(test$p.value)
}

tic("reading inputs:")
###### Input arguments ##############
args = commandArgs(trailingOnly = TRUE)
directory = args[1]
gtf_file = args[2]
is.SE = as.numeric(args[3])
is.10X = as.numeric(args[4])
is.stranded = as.numeric(args[5])
#####################################
toc()

### arguments for debugging ######
#is.SE = 1
#directory = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/test/TSP2_BM_vertebralbody_10X_1_2_5prime_S35_L002/"
#gtf_file =  "/oak/stanford/groups/horence/circularRNApipeline_Cluster/index/grch38_known_genes.gtf"
##################################

###### read in class input file #####################
class_input_file = list.files(directory, pattern = "class_input.tsv")
class_input_file = paste(directory, class_input_file,sep = "")
class_input =  fread(class_input_file, sep = "\t", header = TRUE)
###############################################


# do deduplication for 10X data
if(is.10X == 1){
  class_input = class_input[!duplicated(paste(barcode,UMI,refName_newR1))]
}

setkey(class_input,refName_newR1)
if(is.SE ==0){  # I want to discard those reads that have an unaligned R2 in PE data
  class_input = class_input[!is.na(nmmR2A)]
}

## if data is unstranded (such as smartseq), I want to make reads for A-B and B-A junctions consistent by changing one junction to the other one
if (is.stranded==0){
  tic("unstranded data modification")
  class_input[, numReads:=length(unique(id)), by = refName_newR1]
  class_input[,pos1:=paste(juncPosR1A,juncPosR1B,sep="-"),by=refName_newR1]
  class_input[,pos2:=paste(juncPosR1B,juncPosR1A,sep="-"),by=refName_newR1]
  class_input[pos1%in%class_input$pos2,ambig:=1] # ambig determines that both junctions A-B and B-A are present in the class input
  class_input[(geneR1A_uniq=="unknown") & (geneR1B_uniq=="unknown"),num_unknown_gene:=2] # count the number of unknown gene names in the junction id
  class_input[!((geneR1A_uniq=="unknown") | (geneR1B_uniq=="unknown")),num_unknown_gene:=0]
  class_input[is.na(num_unknown_gene),num_unknown_gene:=1]
  
  ambiguous_juncs = merge(class_input[!is.na(ambig)&!duplicated(refName_newR1),list(refName_newR1,numReads,juncPosR1A,juncPosR1B,num_unknown_gene,pos1)],class_input[!duplicated(refName_newR1),list(refName_newR1,numReads,num_unknown_gene,pos2)],all.x=TRUE,all.y=FALSE,by.x="pos1",by.y="pos2")
  ambiguous_juncs[,swap:=0] # swap determines which id refName_newR1.x or refName_newR1.y should be selected as the final id
  ambiguous_juncs[num_unknown_gene.x > num_unknown_gene.y, swap:=1] # I prefer to select the id that has fewer unknown gene names
  ambiguous_juncs[(num_unknown_gene.x==num_unknown_gene.y) & (numReads.x<numReads.y),swap:=1]
  
  ambiguous_juncs[,max_pos:=max(juncPosR1A,juncPosR1B),by=refName_newR1.x]
  ambiguous_juncs[,min_pos:=min(juncPosR1A,juncPosR1B),by=refName_newR1.x]
  ambiguous_juncs = ambiguous_juncs[!duplicated(paste(max_pos,min_pos))]
  ambiguous_juncs[swap==1,final_refName:=refName_newR1.y] # the selected junction id
  ambiguous_juncs[swap==0,final_refName:=refName_newR1.x]
  ambiguous_juncs[swap==1,old_refName:=refName_newR1.x]  # the discarded junction id
  ambiguous_juncs[swap==0,old_refName:=refName_newR1.y]
  
  class_input = merge(class_input,ambiguous_juncs[,list(final_refName,old_refName)],all.x=TRUE,all.y=FALSE,by.x="refName_newR1",by.y="old_refName")
  class_input[is.na(final_refName),final_refName:=refName_newR1]
  class_input[,refName_newR1:=final_refName]
  class_input[,final_refName:=NULL]
  toc()
}

tic()
class_input[,c("junc_cdf_glm", "junc_cdf_glm_corrected", "junc_cdf_glmnet", "junc_cdf_glmnet_constrained", "junc_cdf_glmnet_corrected", "junc_cdf_glmnet_corrected_constrained"):=NULL]
class_input[,c("p_predicted_glm", "p_predicted_corrected", "p_predicted_glmnet", "p_predicted_glmnet_constrained", "p_predicted_glmnet_corrected", "p_predicted_glmnet_corrected_constrained"):=NULL]
class_input[, numReads :=length(unique(id)), by = refName_newR1]
class_input[, chrR1A:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][1], by = refName_newR1]
class_input[, chrR1B:=strsplit(refName_newR1, split = "[:|]")[[1]][5], by = refName_newR1]
class_input[, juncPosR1A:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][3]), by = refName_newR1]
class_input[, juncPosR1B:=as.integer(strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][6]), by = refName_newR1]
class_input[, gene_strandR1A:=strsplit(refName_newR1, split = "[:|]")[[1]][4], by = refName_newR1]
class_input[, gene_strandR1B:=strsplit(refName_newR1, split = "[:|]")[[1]][8], by = refName_newR1]
class_input[, geneR1A_uniq:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][2], by = refName_newR1]
class_input[, geneR1B_uniq:=strsplit(refName_newR1, split = ":", fixed = TRUE)[[1]][5], by = refName_newR1]
toc()

## add ensembl ids
class_input = add_ensembl(gtf_file,directory,class_input,is.SE)

class_input[fileTypeR1 == "Chimeric",is.STAR_Chim := ""]
class_input[fileTypeR1 == "Aligned",is.STAR_SJ := ""]

## compare class inpout file with STAR chimeric and output files
#class_input = compare_classinput_STARChimOut(directory,is.SE)

### obtain fragment lengths for chimeric reads for computing length adjusted AS ##########
class_input[fileTypeR1 == "Aligned", length_adj_AS_R1:= aScoreR1A/readLenR1]
class_input[fileTypeR1 == "Chimeric", length_adj_AS_R1 := (aScoreR1A + aScoreR1B)/ (MR1A + MR1B + SR1A + SR1B)]
class_input[, length_adj_AS_R1A:= aScoreR1A / (MR1A + SR1A), by = 1:nrow(class_input)]
class_input[, length_adj_AS_R1B:= aScoreR1B / (MR1B + SR1B), by = 1:nrow(class_input)]
###########################################################################################

### obtain the number of smismatches per alignment ##########
class_input[fileTypeR1 == "Aligned", nmmR1:= nmmR1A]
class_input[fileTypeR1 == "Chimeric", nmmR1:= nmmR1A + nmmR1B]
###########################################################################################

class_input[NHR1A>1,multimapping:=1]
class_input[NHR1A==1,multimapping:=0]
class_input[,frac_multimapping:=mean(multimapping),by=refName_newR1]

#### obtain junction overlap #########
class_input[, overlap_R1 := min(MR1A,MR1B), by = 1:nrow(class_input)]
class_input[, max_overlap_R1 := max(MR1A,MR1B), by = 1:nrow(class_input)]
class_input[, median_overlap_R1 := as.integer(round(median(overlap_R1))), by = refName_newR1]
min_overlap_R1 = min(class_input$overlap_R1)
max_overlap_R1 = max(class_input$overlap_R1)
class_input[, sd_overlap:=sqrt(var(overlap_R1)), by = refName_newR1]
######################################

###### compute noisy junction score ########
round.bin = 50
class_input[, binR1A:=round(juncPosR1A/round.bin)*round.bin]
class_input[, binR1B:=round(juncPosR1B/round.bin)*round.bin]
class_input[, njunc_binR1A:=length(unique(refName_newR1)), by = paste(binR1A, chrR1A)]
class_input[, njunc_binR1B:=length(unique(refName_newR1)), by = paste(binR1B, chrR1B)]
class_input[, binR1A:= NULL]
class_input[, binR1B:= NULL]
############################################

## assigning bins and chromosomes to unknown genes
options(scipen = 999)
round.bin = 100000
class_input[geneR1A_uniq=="",geneR1A_uniq:="unknown"]
class_input[geneR1B_uniq=="",geneR1B_uniq:="unknown"]
class_input[geneR1A_uniq=="unknown", binR1A:=round(juncPosR1A/round.bin)*round.bin]
class_input[geneR1B_uniq=="unknown", binR1B:=round(juncPosR1B/round.bin)*round.bin]
class_input[geneR1A_uniq=="unknown",geneR1A_uniq:=paste("unknown",chrR1A,binR1A,sep = "_")]
class_input[geneR1B_uniq=="unknown",geneR1B_uniq:=paste("unknown",chrR1B,binR1B,sep = "_")]
class_input[, binR1A:= NULL]
class_input[, binR1B:= NULL]
############################################

#### get the number of distinct partners for each splice site ###########
setkey(class_input,refName_newR1)
class_input[, junc_pos1_R1:=paste(chrR1A, juncPosR1A, sep= ":"), by = refName_newR1]
class_input[, junc_pos2_R1:=paste(chrR1B, juncPosR1B, sep= ":"), by = refName_newR1]
class_input[, threeprime_partner_number_R1:=length(unique(junc_pos2_R1)), by = junc_pos1_R1]
class_input[, fiveprime_partner_number_R1:=length(unique(junc_pos1_R1)), by = junc_pos2_R1]
class_input[, junc_pos1_R1:= NULL]
class_input[, junc_pos2_R1:= NULL]
###########################################################################


## the same predictors for R2
if (is.SE == 0){
  ### obtain fragment lengths for chimeric reads for computing length adjusted AS ##########
  class_input[fileTypeR2 == "Aligned", length_adj_AS_R2:= aScoreR2A/readLenR2]
  class_input[fileTypeR2 == "Chimeric", length_adj_AS_R2 := (aScoreR2A + aScoreR2B)/ (MR2A + MR2B + SR2A + SR2B)]
  ###########################################################################################
  
  ### obtain the number of mismatches per alignment ##########
  class_input[fileTypeR2 == "Aligned", nmmR2:= nmmR2A]
  class_input[fileTypeR2 == "Chimeric", nmmR2:= nmmR2A + nmmR2B]
  ##########################################################################
  
  #### obtain junction overlap #########
  class_input[, overlap_R2 := min(MR2A,MR2B), by = 1:nrow(class_input)]
  class_input[, max_overlap_R2 := max(MR2A,MR2B), by = 1:nrow(class_input)]
  ######################################
}

class_input[, cur_weight := NULL]
class_input[, train_class := NULL]
####### Assign pos and neg training data for GLM training #######
n.neg = nrow(class_input[genomicAlignmentR1 ==1 & fileTypeR1 == "Aligned"])
n.pos = nrow(class_input[genomicAlignmentR1 ==0 & entropyR1>2 & fileTypeR1 == "Aligned"])
n.neg = min(n.neg,150000)
n.pos = min(n.pos,150000)  # number of positive reads that we want to subsample from the list of all reads
all_neg_reads = which((class_input$genomicAlignmentR1 ==1)  & (class_input$fileTypeR1 == "Aligned"))
all_pos_reads = which((class_input$genomicAlignmentR1 ==0) & (class_input$entropyR1>2) & (class_input$fileTypeR1 == "Aligned"))
class_input[sample(all_neg_reads, n.neg, replace= FALSE), train_class := 0]
class_input[sample(all_pos_reads, n.pos, replace= FALSE), train_class := 1]
#################################################################


#set the training reads and class-wise weights for the reads within the same class
class.weight = min(n.pos, n.neg)

if (n.pos >= n.neg){
  class_input[train_class == 0, cur_weight := 1]
  class_input[train_class == 1, cur_weight := n.neg / n.pos]
} else {
  class_input[train_class == 0, cur_weight := n.pos/n.neg]
  class_input[train_class == 1, cur_weight := 1]
}
####################################


####################################################
######### GLMnet model (constrained)  ##############
####################################################
tic("GLMnet constrained")
print("GLMnet constrained")

if (is.SE == 0){
  regression_formula = as.formula("train_class ~ overlap_R1 * max_overlap_R1 + NHR1A + nmmR1 + MR1A:SR1A + MR1B:SR1B + length_adj_AS_R1 + nmmR2 + length_adj_AS_R2 + NHR2A + entropyR1*entropyR2 + location_compatible + read_strand_compatible")
} else {
  regression_formula = as.formula("train_class ~ overlap_R1 * max_overlap_R1 + NHR1A + nmmR1 + MR1A:SR1A +  MR1B:SR1B + entropyR1 + length_adj_AS_R1 + entropyR1:NHR1A + entropyR1:length_adj_AS_R1")
}

x_glmnet = model.matrix(regression_formula, class_input[!is.na(train_class)])
if (is.SE==0){
  glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(class_input[!is.na(train_class)]$train_class), family =c("binomial"), class_input[!is.na(train_class)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits=c(rep(Inf,3),0,0, rep(Inf,12)), lower.limits=c(-Inf,0,0, rep(-Inf,2),0, rep(-Inf,3),0,0,0,0, rep(-Inf,4)) )
}else{
  glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(class_input[!is.na(train_class)]$train_class), family =c("binomial"), class_input[!is.na(train_class)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits=c(rep(Inf,3),0,0, rep(Inf,7)), lower.limits=c(-Inf,0,0, rep(-Inf,2),0,0, rep(-Inf,4),0) )
}
print("done with fitting GLMnet constrained")
toc()
print(coef(glmnet_model_constrained, s = "lambda.1se"))
print(glmnet_model_constrained)

# predict for all read alignments in the class input file
class_input_glmnet = model.matrix(update(regression_formula, refName_newR1 ~ .), class_input)

if(nrow(class_input_glmnet)>20000000){ # if class input is too large, we apply the predict function sequentially on the smaller chunks of the data frame
  num_iter = ceiling(nrow(class_input_glmnet)/20000000)
  a = c()
  for (counter in 1:num_iter){
    start = (counter-1)*20000000 + 1
    end = min(counter*20000000, nrow(class_input_glmnet))
    a[start:end] = predict(glmnet_model_constrained, newx = class_input_glmnet[start:end,], type = "response", s = "lambda.1se", se.fit = TRUE)
  }
  class_input$glmnet_per_read_prob_constrained = a
} else{
  class_input$glmnet_per_read_prob_constrained = predict(glmnet_model_constrained, newx = class_input_glmnet, type = "response", s = "lambda.1se", se.fit = TRUE)
}


# compute the fitted classification error based on training data
compute_class_error(class_input[!is.na(train_class)]$train_class, class_input[!is.na(train_class)]$glmnet_per_read_prob_constrained)

# compute aggregated score for each junction
class_input[, p_predicted_glmnet_constrained:= 1/( exp(sum(log( (1 - glmnet_per_read_prob_constrained)/glmnet_per_read_prob_constrained ))) + 1), by = refName_newR1]


# compute the junc_cdf scores
class_input = compute_junc_cdf(class_input , "p_predicted_glmnet_constrained", "glmnet_per_read_prob_constrained", "junc_cdf_glmnet_constrained")
print("done with GLMnet constrained")


if (is.SE==0){
  tic("GLMnet corrected constrained")
  class_input[, glmnet_per_read_prob_corrected_constrained := glmnet_per_read_prob_constrained]
  class_input[(location_compatible==0 | read_strand_compatible==0), glmnet_per_read_prob_corrected_constrained:=glmnet_per_read_prob_constrained/(1 + glmnet_per_read_prob_constrained)]
  class_input[, p_predicted_glmnet_corrected_constrained := 1/( exp(sum(log( (1 - glmnet_per_read_prob_corrected_constrained)/glmnet_per_read_prob_corrected_constrained ))) + 1), by = refName_newR1]
  class_input = compute_junc_cdf(class_input , "p_predicted_glmnet_corrected_constrained", "glmnet_per_read_prob_corrected_constrained", "junc_cdf_glmnet_corrected_constrained")
  print("done with GLMnet corrected contrained")
  toc()  
}


######################################################
######################################################
#### now we do the two-step GLM for chimeric reads ###
######################################################
######################################################
if (nrow(class_input[fileTypeR1=="Chimeric"])>0){
  class_input[, train_class := NULL]
  class_input[, cur_weight := NULL]
  p_predicted_quantile = quantile(class_input[!(duplicated(refName_newR1)) & (fileTypeR1 == "Chimeric")]$p_predicted_glmnet_constrained, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) # the percentiles for the scores based on linear GLM
  p_predicted_neg_cutoff = p_predicted_quantile[[2]]
  p_predicted_pos_cutoff = p_predicted_quantile[[8]]
  
  ####### Assign pos and neg training data for GLM training #######
  n.neg = nrow(class_input[fileTypeR1 == "Chimeric"][(p_predicted_glmnet_constrained <= p_predicted_neg_cutoff) | (refName_newR1%like%"chrM")])
  n.pos = nrow(class_input[fileTypeR1 == "Chimeric"][(p_predicted_glmnet_constrained >= p_predicted_pos_cutoff) & (!refName_newR1%like%"chrM")])
  n.neg = min(n.neg,150000)
  n.pos = min(n.pos,150000)  # number of positive reads that we want to subsample from the list of all reads
  all_neg_reads = which((class_input$fileTypeR1 == "Chimeric") & ((class_input$p_predicted_glmnet_constrained <= p_predicted_neg_cutoff) | class_input$refName_newR1%like%"chrM"))
  all_pos_reads = which((class_input$fileTypeR1 == "Chimeric") & ((class_input$p_predicted_glmnet_constrained >= p_predicted_pos_cutoff) & !class_input$refName_newR1%like%"chrM"))
  class_input[sample(all_neg_reads, n.neg, replace= FALSE), train_class := 0]
  class_input[sample(all_pos_reads, n.pos, replace= FALSE), train_class := 1]
  #################################################################
  
  class.weight = min(n.pos, n.neg)
  if (n.pos >= n.neg){
    class_input[train_class == 0, cur_weight := 1]
    class_input[train_class == 1, cur_weight := n.neg / n.pos]
  } else {
    class_input[train_class == 0, cur_weight := n.pos/n.neg]
    class_input[train_class == 1, cur_weight := 1]
  }
  
  
  ################################################################
  ######## Building the GLMnet for chimeric junctions ############
  ################################################################
  if (is.SE == 0){
    regression_formula = as.formula("train_class ~ overlap_R1 * max_overlap_R1  + nmmR1 + length_adj_AS_R1A + length_adj_AS_R1B + nmmR2 + entropyR1*entropyR2 + length_adj_AS_R2")
  } else{
    regression_formula = as.formula("train_class ~ overlap_R1 * max_overlap_R1  + nmmR1 + length_adj_AS_R1A + length_adj_AS_R1B + entropyR1")
  }
  
  tic("Two-step GLMnet constrained")
  print("Two-step GLMnet constrained")
  x_glmnet = model.matrix(regression_formula, class_input[!is.na(train_class)])
  if (is.SE == 0){
    glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(class_input[!is.na(train_class)]$train_class), family =c("binomial"), class_input[!is.na(train_class)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits = c(rep(Inf,3),0, rep(Inf,8)), lower.limits = c(-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,-Inf,-Inf))
  }else{
    glmnet_model_constrained = cv.glmnet(x_glmnet, as.factor(class_input[!is.na(train_class)]$train_class), family =c("binomial"), class_input[!is.na(train_class)]$cur_weight, intercept = FALSE, alpha = 1, nlambda = 50, nfolds = 5, upper.limits = c(rep(Inf,3),0, rep(Inf,4)), lower.limits = c(-Inf,0,0,-Inf,0,0,0,-Inf))
  }
  print("done with fitting Two-step GLMnet constrained")
  print(coef(glmnet_model_constrained, s = "lambda.1se"))
  toc()
  
  # predict for all chimeric alignments in the class input file
  class_input_glmnet_constrained = model.matrix(update(regression_formula, refName_newR1 ~ .), class_input[fileTypeR1 == "Chimeric"])
  class_input[,glmnet_twostep_per_read_prob_constrained:=0]
  a = predict(glmnet_model_constrained, newx = class_input_glmnet_constrained, type = "response", s = "lambda.1se", se.fit = TRUE)
  class_input[fileTypeR1 == "Chimeric",glmnet_twostep_per_read_prob_constrained:=a]
  
  # compute the fitted classification error based on training data
  compute_class_error(class_input[!is.na(train_class)]$train_class, class_input[!is.na(train_class)]$glmnet_twostep_per_read_prob_constrained)
  
  # compute aggregated score for each junction
  class_input[fileTypeR1 == "Chimeric", p_predicted_glmnet_twostep_constrained:= 1/( exp(sum(log( (1 - glmnet_twostep_per_read_prob_constrained)/glmnet_twostep_per_read_prob_constrained ))) + 1), by = refName_newR1]
  ######################################
  ######################################
}

setkey(class_input,refName_newR1)
class_input[, frac_genomic_reads :=mean(genomicAlignmentR1), by = refName_newR1]

class_input[, ave_AT_run_R1:=mean(AT_run_R1), by = refName_newR1]
class_input[, ave_GC_run_R1:=mean(GC_run_R1), by = refName_newR1]
class_input[, ave_max_run_R1:=mean(max_run_R1), by = refName_newR1]

class_input[, ave_entropyR1:=mean(entropyR1), by = refName_newR1]
class_input[, min_entropyR1:=min(entropyR1), by = refName_newR1]

if (is.SE == 0){
  class_input[, frac_anomaly:=0]
  class_input[(location_compatible==0 | read_strand_compatible==0), frac_anomaly:=.N/numReads, by = refName_newR1] # the fraction of anomalous reads for each junction
  class_input[, ave_AT_run_R2:=mean(AT_run_R2), by = refName_newR1]
  class_input[, ave_GC_run_R2:=mean(GC_run_R2), by = refName_newR1]
  class_input[, ave_max_run_R2:=mean(max_run_R2), by = refName_newR1]
  class_input[, ave_entropyR2:=mean(entropyR2), by = refName_newR1]
  class_input[, min_entropyR2:=min(entropyR2), by = refName_newR1]
}

###############################################################################
######### p-value for junction median overlap #################################
### compute p-value for how close is the median overlap to what we expect based on the read length
iter=5000
tic("junc_median_p_val")
for (num_reads in 1:15){
  rnd_overlaps = matrix(0, iter, num_reads)
  rnd_overlaps = apply(rnd_overlaps,1, function(x) sample(min_overlap_R1:max_overlap_R1, num_reads, replace = TRUE))
  rnd_overlaps = t(rnd_overlaps)
  if(num_reads == 1){
    rnd_overlaps = t(rnd_overlaps)  # for num_reads=1 I need to transepose twice since first I have a vector
  }
  null_dist_medians = apply(rnd_overlaps,1, function(x) median(x))
  class_input[numReads == num_reads, p_val_median_overlap_R1:=length(which(null_dist_medians > median_overlap_R1))/iter, by = median_overlap_R1]
}
class_input[numReads > 15, p_val_median_overlap_R1:=pnorm(median_overlap_R1, mean = (min_overlap_R1+max_overlap_R1)/2, sd = sqrt( (max_overlap_R1-min_overlap_R1)^2 /12 /numReads), lower.tail = FALSE), by = refName_newR1]
toc()
#####################################



col_names_to_keep_in_junc_pred_file = c("refName_newR1","frac_genomic_reads","frac_multimapping","numReads","njunc_binR1B","njunc_binR1A","median_overlap_R1","threeprime_partner_number_R1","fiveprime_partner_number_R1","is.STAR_Chim","is.STAR_SJ","is.STAR_Fusion","is.True_R1","geneR1A_expression_stranded","geneR1A_expression_unstranded","geneR1B_expression_stranded","geneR1B_expression_unstranded","geneR1A_RPKM_stranded","geneR1A_RPKM_unstranded","geneR1B_RPKM_stranded","geneR1B_RPKM_unstranded","geneR1B_ensembl","geneR1A_ensembl","geneR1B_uniq","geneR1A_uniq","intron_motif","is.TRUE_fusion","p_predicted_glmnet_constrained","p_predicted_glmnet_corrected_constrained","p_predicted_glmnet_twostep_constrained","junc_cdf_glmnet_constrained","junc_cdf_glmnet_corrected_constrained","junc_cdf_glmnet_twostep","ave_max_junc_14mer","ave_min_junc_14mer","frac_anomaly","ave_AT_run_R1","ave_GC_run_R1","ave_max_run_R1","ave_AT_run_R2","ave_GC_run_R2","ave_entropyR1","ave_entropyR2","min_entropyR1","min_entropyR2","ave_max_run_R2","sd_overlap","p_val_median_overlap_R1","chrR1A","chrR1B","juncPosR1A","juncPosR1B","gene_strandR1A","gene_strandR1B","fileTypeR1","read_strandR1A","read_strandR1B")
GLM_output = unique(class_input[, colnames(class_input)%in%col_names_to_keep_in_junc_pred_file, with = FALSE])
GLM_output = GLM_output[!(duplicated(refName_newR1))]

##############################################################################################################
#### compute emp.p values based upon junc_cdf and using junctions with >10% genomic reads ####################
null_dist = GLM_output[is.na(is.STAR_Chim) & frac_genomic_reads > 0.1]$junc_cdf_glmnet_constrained
GLM_output[, emp.p_glmnet_constrained:=length(which(null_dist>junc_cdf_glmnet_constrained))/length(null_dist), by = junc_cdf_glmnet_constrained]

if (is.SE == 0){
  null_dist = GLM_output[is.na(is.STAR_Chim) & frac_genomic_reads>0.1]$junc_cdf_glmnet_corrected_constrained
  GLM_output[, emp.p_glmnet_corrected_constrained:=length(which(null_dist>junc_cdf_glmnet_corrected_constrained))/length(null_dist), by = junc_cdf_glmnet_corrected_constrained]
}
##############################################################################################################
##############################################################################################################


##########################################################
##### add junction sequence to the glm report file  ######
tic("add_junc_seq")
class_input_extract = class_input[, list(refName_newR1, seqR1,  flagR1A, flagR1B,cigarR1A, cigarR1B)]
setkey(class_input_extract,cigarR1A)
class_input_extract[, readoverhang1_length:=sum(explodeCigarOpLengths(cigarR1A, ops=c("I", "S","M"))[[1]]), by = cigarR1A]
setkey(class_input_extract,cigarR1B)
class_input_extract[, readoverhang2_length:=sum(explodeCigarOpLengths(cigarR1B, ops=c("I", "S","M"))[[1]]), by = cigarR1B]
class_input_extract[, overlap:=min(readoverhang1_length, readoverhang2_length), by = 1:nrow(class_input_extract)]   # I do these steps as I want to get the junction sequence based on the most balanced read alignment for the junction (the one that has the highest overlap) 
class_input_extract[, max_junc_level_overlap:=max(overlap), by = refName_newR1]  # the maximum overlap across all reads aligned to the junction 
class_input_extract = class_input_extract[overlap==max_junc_level_overlap]
class_input_extract[, c("overlap","max_junc_level_overlap"):= NULL]
class_input_extract = class_input_extract[!duplicated(refName_newR1)]
GLM_output = merge(GLM_output, class_input_extract, by.x = "refName_newR1", by.y = "refName_newR1", all.x = TRUE, all.y = FALSE)
toc()
##########################################################
##########################################################


## removing redundant columns for GLM script
class_input[,c("cur_weight","train_class","sum_log_per_read_prob","log_per_read_prob","junc_cdf1_glmnet_twostep","refName_readStrandR1","refName_readStrandR2","gene_strandR1A_new","gene_strandR1B_new"):=NULL]
GLM_output[, c("cigarR1A","cigarR1B"):= NULL]

GLM_output = postprocessing_model(GLM_output,class_input,is.SE,is.10X)

write.table(GLM_output, paste(directory,"GLM_output.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(class_input, paste(directory,"class_input_tmp.tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
system( paste("mv ",directory,"class_input_tmp.tsv ",directory,"class_input.tsv ",sep=""))

# use the following function to call final splice junctions 
sicilian_splicing_called_junctions = postprocessing_hardthreshold(GLM_output,class_input,is.10X,is.SE)


# use the following function to find inserted/missing domains for each called junction (this function will be run only when ucsc domain annotation files have been provided)
if ((length(args)==6) | (length(args)==8)){
  ucsc_domain_file = args[6]
  sicilian_splicing_called_junctions = domain_annotation(ucsc_domain_file,sicilian_splicing_called_junctions)
}


write.table(sicilian_splicing_called_junctions, paste(directory,"sicilian_called_splice_juncs.tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")

# use this for annotating the exon boudnaries in the called junctions (this function will be run only when exon and splice annotation pickle files have been provided)
if ((length(args)==7) | (length(args)==8)){
  exon_pickle = args[length(args)-1]
  splice_pickle = args[length(args)]
  script_directory = getwd()
  system(paste("python3 ",script_directory,"/scripts/ann_splices.py -i ",directory,"sicilian_called_splice_juncs.tsv", " -o ",directory,"sicilian_called_splice_juncs.tsv"," -e ",exon_pickle, " -s ", splice_pickle,sep = ""))
}
toc()
