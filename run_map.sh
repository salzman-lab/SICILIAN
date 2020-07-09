#!/bin/bash
#
#SBATCH --job-name=map_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003
#SBATCH --output=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/map_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.out
#SBATCH --error=/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/log_files/map_lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003.%j.err
#SBATCH --time=24:00:00
#SBATCH --account=horence
#SBATCH --partition=nih_s10
#SBATCH --nodes=1
#SBATCH --mem=60Gb
date
mkdir -p /oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003
STAR --version
/oak/stanford/groups/horence/Roozbeh/single_cell_project/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 4 --genomeDir /oak/stanford/groups/horence/Roozbeh/single_cell_project/STAR-2.7.1a/hg38_index_2.7.1a --readFilesIn /oak/stanford/groups/krasnow/ktrav/COVID/data10x/sequencing_runs/200603_A00111_0499_BH7WTJDSXY/fastqs/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003_extracted_R2_001.fastq.gz --readFilesCommand zcat --twopassMode Basic --alignIntronMax 1000000 --alignIntronMin 20 --outFileNamePrefix /oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Krasnow_COVID_10x_pilot3/lungSlice_Pilot3_72h_SARS-CoV-2_Sample_2_S8_L003/2 --outSAMtype BAM Unsorted --chimSegmentMin 10 --outSAMattributes All --chimOutType WithinBAM HardClip Junctions --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignMatesGapMax 1000000 --scoreGapNoncan -4 --scoreGapATAC -4 --chimScoreJunctionNonGTAG 0 --limitOutSJcollapsed 5000000 --limitIObufferSize 250000000 --chimJunctionOverhangMin 10 --alignSJstitchMismatchNmax -1 -1 -1 -1 --chimSegmentReadGapMax 0--quantMode GeneCounts --limitSjdbInsertNsj 2600000 --sjdbGTFfile /oak/stanford/groups/krasnow/ktrav/COVID/data10x/sequencing_runs/200528_A00111_0493_BHLJ7KDRXX/gencode-vH29.SARS-CoV-2_WA1.gtf --outReadsUnmapped Fastx 


date
