#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(plotROC)

## this script is used for generating the ROC curves (Figure 1E in the SICILIAN paper, the first ) based on the simulated benchmarking datrasets,
## when splice junctions are called based on the SICILIAN or read count criterion

working_directory = getwd()

########################################################################################################
######## The ROC curve by SICILIAN for HISAT benchmarking data (perfect dataset)  ######################
########################################################################################################
glm_output = fread(paste(working_directory,"/benchmarking_files/HISAT/HISAT_perfect/GLM_output_benchmark_OL1.txt",sep=""),sep="\t",header=TRUE) #output of the glm script

glm_output[,SICILIAN:=1-emp.p_glmnet_corrected_constrained]
setnames(glm_output,"numReads","Read_count")
data_for_ROC <- melt_roc(glm_output, "is.True_R1", c("SICILIAN", "Read_count"))

#### below we compute the AUC values based on the SICILIAN and read count criteria  #######
ROC_plot_SICILIAN = ggplot()+geom_roc(aes(d = is.True_R1, m = SICILIAN, color = "SICILIAN"), glm_output)
AUC_SICILIAN = round(calc_auc(ROC_plot_SICILIAN)$AUC, 3) # AUC value based on SICILIAN
ROC_plot_readcount = ggplot()+geom_roc(aes(d = is.True_R1, m = Read_count, color = "Read count"), glm_output)
AUC_read_count = round(calc_auc(ROC_plot_readcount)$AUC, 3) # AUC values based on read count criterion
###########################################################################################

## plot the ROC curves for both SICILIAN and read count in one figure along with their AUC values
ROC_plot =  ggplot(data_for_ROC , aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
ROC_plot + annotate(geom = "text", x = .75, y = .5, label = paste("AUC_SICILIAN =", AUC_SICILIAN)) + annotate(geom = "text", x = .75, y = .25, label = paste("AUC_read_count =", AUC_read_count)) + coord_fixed(ratio = 0.6)+theme(axis.text=element_text(size=12), 
                             axis.title=element_text(size = 14,face = "bold"),legend.text = element_text(size = 16))
########################################################################################################
########################################################################################################
########################################################################################################


########################################################################################################
######## The ROC curve by SICILIAN for HISAT benchmarking data (mismatch dataset) ######################
########################################################################################################
glm_output = fread(paste(working_directory,"/benchmarking_files/HISAT/HISAT_mismatch/GLM_output_benchmark_OL1.txt",sep=""),sep="\t",header=TRUE) #output of the glm script

glm_output[,SICILIAN:=1-emp.p_glmnet_corrected_constrained]
setnames(glm_output,"numReads","Read_count")
data_for_ROC <- melt_roc(glm_output, "is.True_R1", c("SICILIAN", "Read_count"))

#### below we compute the AUC values based on the SICILIAN and read count criteria  #######
ROC_plot_SICILIAN = ggplot()+geom_roc(aes(d = is.True_R1, m = SICILIAN, color = "SICILIAN"), glm_output)
AUC_SICILIAN = round(calc_auc(ROC_plot_SICILIAN)$AUC, 3) # AUC value based on SICILIAN
ROC_plot_readcount = ggplot()+geom_roc(aes(d = is.True_R1, m = Read_count, color = "Read count"), glm_output)
AUC_read_count = round(calc_auc(ROC_plot_readcount)$AUC, 3) # AUC values based on read count criterion
###########################################################################################

## plot the ROC curves for both SICILIAN and read count in one figure along with their AUC values
ROC_plot =  ggplot(data_for_ROC , aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
ROC_plot + annotate(geom = "text", x = .75, y = .5, label = paste("AUC_SICILIAN =", AUC_SICILIAN)) + annotate(geom = "text", x = .75, y = .25, label = paste("AUC_read_count =", AUC_read_count)) + coord_fixed(ratio = 0.6)+theme(axis.text=element_text(size=12),
                              axis.title=element_text(size = 14,face = "bold"),legend.text = element_text(size = 16))
########################################################################################################
########################################################################################################
########################################################################################################


################################################################################
######## The ROC curve by SICILIAN for Engstrom data sim1 ######################
################################################################################
glm_output = fread(paste(working_directory,"/benchmarking_files/Engstrom/sim1/GLM_output_benchmark_OL1.txt",sep=""),sep="\t",header=TRUE) #output of the glm script

glm_output[,SICILIAN:=1-emp.p_glmnet_corrected_constrained]
setnames(glm_output,"numReads","Read_count")
data_for_ROC <- melt_roc(glm_output, "is.True_R1", c("SICILIAN", "Read_count"))

#### below we compute the AUC values based on the SICILIAN and read count criteria  #######
ROC_plot_SICILIAN = ggplot()+geom_roc(aes(d = is.True_R1, m = SICILIAN, color = "SICILIAN"), glm_output)
AUC_SICILIAN = round(calc_auc(ROC_plot_SICILIAN)$AUC, 3) # AUC value based on SICILIAN
ROC_plot_readcount = ggplot()+geom_roc(aes(d = is.True_R1, m = Read_count, color = "Read count"), glm_output)
AUC_read_count = round(calc_auc(ROC_plot_readcount)$AUC, 3) # AUC values based on read count criterion
###########################################################################################

## plot the ROC curves for both SICILIAN and read count in one figure along with their AUC values
ROC_plot =  ggplot(data_for_ROC , aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
ROC_plot + annotate(geom = "text", x = .75, y = .5, label = paste("AUC_SICILIAN =", AUC_SICILIAN)) + annotate(geom = "text", x = .75, y = .25, label = paste("AUC_read_count =", AUC_read_count)) + coord_fixed(ratio = 0.6)+theme(axis.text=element_text(size=12), 
                                    axis.title=element_text(size = 14,face = "bold"),legend.text = element_text(size = 16))
##################################################################################
##################################################################################
##################################################################################


################################################################################
######## The ROC curve by SICILIAN for Engstrom data sim2 ######################
################################################################################
glm_output = fread(paste(working_directory,"/benchmarking_files/Engstrom/sim2/GLM_output_benchmark_OL1.txt",sep=""),sep="\t",header=TRUE) #output of the glm script

glm_output[,SICILIAN:=1-emp.p_glmnet_corrected_constrained]
setnames(glm_output,"numReads","Read_count")
data_for_ROC <- melt_roc(glm_output, "is.True_R1", c("SICILIAN", "Read_count"))

#### below we compute the AUC values based on the SICILIAN and read count criteria  #######
ROC_plot_SICILIAN = ggplot()+geom_roc(aes(d = is.True_R1, m = SICILIAN, color = "SICILIAN"), glm_output)
AUC_SICILIAN = round(calc_auc(ROC_plot_SICILIAN)$AUC, 3) # AUC value based on SICILIAN
ROC_plot_readcount = ggplot()+geom_roc(aes(d = is.True_R1, m = Read_count, color = "Read count"), glm_output)
AUC_read_count = round(calc_auc(ROC_plot_readcount)$AUC, 3) # AUC values based on read count criterion
###########################################################################################

## plot the ROC curves for both SICILIAN and read count in one figure along with their AUC values
ROC_plot =  ggplot(data_for_ROC , aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
ROC_plot + annotate(geom = "text", x = .75, y = .5, label = paste("AUC_SICILIAN =", AUC_SICILIAN)) + annotate(geom = "text", x = .75, y = .25, label = paste("AUC_read_count =", AUC_read_count)) + coord_fixed(ratio = 0.6)+theme(axis.text=element_text(size=12), 
                                  axis.title=element_text(size = 14,face = "bold"),legend.text = element_text(size = 16))
##################################################################################
##################################################################################
##################################################################################