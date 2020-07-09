#!/bin/bash
#
#SBATCH --job-name=class_input_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003
#SBATCH --output=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/class_input_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.out
#SBATCH --error=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/class_input_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.err
#SBATCH --time=48:00:00
#SBATCH --account=horence
#SBATCH --partition=nih_s10
#SBATCH --nodes=1
#SBATCH --mem=300Gb
#SBATCH --dependency=afterok:15949994
#SBATCH --kill-on-invalid-dep=yes
date
python3 scripts/light_class_input.py --outpath /oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/ --gtf /oak/stanford/groups/krasnow/ktrav/COVID/data10x/sequencing_runs/200528_A00111_0493_BHLJ7KDRXX/gencode-vH29.SARS-CoV-2_WA1.gtf --annotator /oak/stanford/groups/horence/Roozbeh/single_cell_project/scripts/STAR_wrapper/annotators/gencode-vH29.SARS-CoV-2_WA1.pkl --bams /oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/2Aligned.out.bam --UMI_bar --stranded_library 
date
