#!/bin/bash
#
#SBATCH --job-name=GLM_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003
#SBATCH --output=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/GLM_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.out
#SBATCH --error=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/GLM_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.err
#SBATCH --time=12:00:00
#SBATCH --account=horence
#SBATCH --partition=nih_s10
#SBATCH --nodes=1
#SBATCH --mem=200Gb
#SBATCH --dependency=afterok:15941426:15941427
#SBATCH --kill-on-invalid-dep=yes
date
Rscript scripts/GLM_script_light.R /oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/ /oak/stanford/groups/krasnow/ktrav/COVID/data10x/sequencing_runs/200528_A00111_0493_BHLJ7KDRXX/gencode-vH29.SARS-CoV-2_WA1.gtf  1 /oak/stanford/groups/horence/Roozbeh/single_cell_project/utility_files/ucscGenePfam.txt /oak/stanford/groups/horence/Roozbeh/single_cell_project/scripts/STAR_wrapper/annotators/gencode-vH29.SARS-CoV-2_WA1_exon_bounds.pkl /oak/stanford/groups/horence/Roozbeh/single_cell_project/scripts/STAR_wrapper/annotators/hg38_RefSeq_splices.pkl 
date
