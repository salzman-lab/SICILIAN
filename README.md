# SICILIAN
![Image of SICILIAN](https://github.com/salzmanlab/SICILIAN/blob/master/SICILIAN.png)

SICILIAN is a statistical method for identifying RNA splice junctions using alignments reported from a spliced aligner. SICILIAN is currently implemented for the STAR aligner, and will be adapted to more spliced aligner in the newar future.

## How to run SICILIAN
There are two options for running SICILIAN
1. using the ready-to-run cloud-based online tool with CWL implementation without the need for preinstallations on the cancer genomics cloud platform.
2. Downloading SICILIAN scripts to a local high performance computing (HPC) cluster and running SICILIAN after installing required libraries. 

## Online cloud-based tool for SICILIAN
The cloud-based computational workflow for SICILIAN has been implemented on the Seven Bridges Cancer Genomics Cloud platform sponsored by the National Cancer Institute. The workflow is fully dockerized and has a user-friendly interface, which facilitates running SICILIAN even for users with little bioinformatics expertise. User only needs to upload his/her own RNA-Seq fastq files or use RNA-Seq files publicly-available on CGC (i.e., TCGA, CCLE, TARGET, ...) to run SICILIAN and selects correct annotation and index files to run SICILIAN. The online tool can be accessed via: https://cgc.sbgenomics.com/u/anaDsbg/sicilian-dev/ 

## Running SICILIAN scripts on a local cluster

- Download the latest version of SICILIAN codes by cloning its github repository 
```
git clone https://github.com/salzmanlab/SICILIAN.git
```
- Intall needed packages for `R` and `Python`
- Build annotation pickle files and STAR index files for a genome assembly and annotation using `create_annotator.py` script
- Set the input variables in `SICILIAN.py` which is the main script that submits all necessary jobs for running SICILIAN on input RNA-Seq data. 
- For running SICILIAN on 10x scRNA-Seq data, it needs to first demultiplex 10x fastq file based on cell barcode and UMI information stored in R1. For 10x analysis, SICILIAN executes `whitelise` and `extract` commands from `UMI_tools` software. Therefore, `UMI_tools` should be preinstalled on the local cluster for running SICILIAN on 10x samples.  
### Software requirements 
- SICILIAN has been developed using `Python 3.6.1` and the following libraries are needed to be installed
  - argparse
  - collections
  - datetime
  - math
  - numpy
  - os
  - pandas
  - pyarrow
  - pickle
  - pysam
  - re
  - time
  - sys

- SICILIAN has been developed with `R 3.6.1` and needs the following packages are needed to be installed in `R`:
  - data.table
  - glmnet
  - tictoc
  - dplyr
  - stringr
  - GenomicAlignments
  (R scripts in SICILIAN automatically check to see if an R library is already installed and then install those that are needed. So no need for manual preinstallation!) 
  
### Building annotator and index files
SICLIAN uses STAR as the aligner and therefore STAR index files need to be input to SICILIAN. STAR index files can be built using the following command:
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir Genome_data/star --genomeFastaFiles fasta_file.fa --sjdbGTFfile gtf_file.gtf
```
SICILIAN also needs annotator pickle files for pulling gene names and adding them to the detected junction based on their coordinates. Annotator pickle file need to be made once for an annotation GTF file using `create_annotator.py` script in the SICILIAN directory. The annotator pickle file can be built with the following command:
```
python3 create_annotator.py -g gtf_file.gtf -a annotation_name
```
`annotation_name` can be set to any arbitrary name but we recommended that it contains the name and version of the annotation (i.e., `hg38_gencode_v33`).    
After running the above command,`create_annotator.py` will create 3 different pickle files: `annotation_name.pkl`, `annotation_name_exon_bounds.pkl`, and `annotation_name_splices.pkl`. 
- `annotation_name.pkl`: is a required input for SICILIAN and is used to add gene names to junction ids
- `annotation_name_exon_bounds.pkl`: is an optional input for SICILIAN and is used to determine whether or not the splice sites in a junction are annotated exon boundaries
- `annotation_name_splices.pkl`: is an optional input for SICILIAN and is used to determine whether or not the splice site is annotated in the annotation file
### Input parameters for SICILIAN:
Update the input parameters in `SICILIAN.py` script with information about your sample, genome assembly and annotations, and STAR alignment.
* `data_path`: specifies path to the directory that contains the fastq files for the input RNA-Seq data.
* `names`: specifies the name of the fastq file for the input RNA-Seq data (without suffix)  
* `r_ends`: list of unique endings for the file names of R1 and R2 fastq files. Example: `r_ends = ["_1.fastq.gz", "_2.fastq.gz"]`
* `out_dir`: specifies path to the directory that will contain the folder specified by `run_name` for SICILIAN output files.
* `run_name`: folder name for the SICILIAN output files.
* `star_path`: the path to the STAR executable file
* `star_ref_path`: the path to t:q
he STAR index files 
* `gtf_file`: the path to the GTF file used as the reference annotation file for the genome assembly.
* `annotator_file`: the path to the `annotation_name.pkl` file
* `exon_pickle_file`: the path to the `annotation_name_exon_bounds.pkl` file (this is an OPTIONAL input)
* `splice_pickle_file`: the path to the `annotation_name_splices.pkl` file (this is an OPTIONAL input)
* `domain_file`: the path to the reference file for annotated protein domains downloaded from UCSC used for finding missing and inserted protein domains in the splice junction (this is an OPTIONAL input)    
* `single`: set equal to `True` if the data is single-end, and `False` if it is paired-end. Note that currently if `single = True` it is assumed that the single read to be aligned is in the second fastq file (because of the tendancy of SICILIAN for droplet (10x) single-cell protocols in which `R1` contains the cell barcode and UMI information and R2 contains the actual cDNA information). This also causes the files to be demultiplexed to create a new fastq file before they're mapped.
* `tenX`: set equal to `True` if the input RNA-Seq data is 10x and `False` otherwise.
* `stranded_library`: set equal to `True` if input RNA-Seq data is based on a stranded library and `False` otherwise. (for stranded libraries such as 10x, `stranded_library` should be set to `True`). When `stranded_library` is set to `True`, strand orientations from the alignment bam file will be used as the strand orientation of the junction. For unstranded libraries, SICILIAN uses gene strand information from the GTF file as the read strand is ambiguous.  
* `bc_pattern`: this parameter is needed only for 10x data and determines the barcode/UMI pattern in R1. For V3 chemistry in which UMI has 12 bps, `bc_pattern` should be set to `"C"*16 + "N"*12` and for 10x data based on V2 chemistry it should be set to `"C"*16 + "N"*12`. `bc_pattern` is needed for `UMI_tools` steps before STAR alignment on input 10x data.
  
### Choosing STAR parameters
STAR alignment parameters can be adjusted in the `STAR_map` function in `SICILIAN.py`. By default, SICILIAN runs STAR with default parameters. 

### Choosing which scripts to run
These parameters let you decide which steps of the SICILIAN pipeline you want to run (for example, if you have run STAR alignment once on a data set but SICILIAN has failed to produce output files, you do not need to run STAR alignment once again and can just set `run_map` to `False` to skip the alignment step when rerunning SICILIAN):
* `run_whitelist`: Set equal to `True` if you want to run UMI-tools whitelist script to extract cell barcodes and identify the most likely true cell barcodes (will be run only for 10X)
* `run_extract`: Set equal to `True` if you want to run UMI-tools extract script which removes UMIs from fastq reads and append them to read name (will be ron only for 10x)
* `run_map`: Set equal to `True` if you want to run the STAR alignment, and `False` otherwise
* `run_class`: Set equal to `True` if you want to run the class input job, and `False` otherwise
* `run_GLM`: Set equal to `True` if you want to run the GLM step and assign statistical scores to each junction in the class input file. The output of this step is a file named `GLM_output.txt`. 

After assigning these variables, run `python3 SICILIAN.py` on the command line to submit the SICILIAN jobs for the input data.

## Description of output files

The output files by SICILIAN  will be written in the output folder specified by `out_dir/run_name` and will contain both the output files generated by STAR and also output files generated by STAR postsprocessing modeling.

### STAR output files
We currently run STAR with `--outSAMtype BAM Unsorted` and `--chimOutType WithinBAM SoftClip Junctions`; therefore, STAR will produce a single BAM file that will contain both Aligned and Chiemric alignments. In addition to the alignment BAM file and STAR log files, there will be `Chimeric.out.junction`, `SJ.out.tab`, and `Unmapped.out.mate` files, containing spliced junctions, summary of all chimeric alignments, and unmapped reads, respectively. For paired-end data as SICILIAN runs STAR independently for `R1` and `R2` reads, STAR output file names for `R1` begin with 1 (i.e., `1Aligned.out.bam`) and those for `R2` begin with 2 (i.e., `1Aligned.out.bam`). For single-end reads, all STAR output files withh begin with 2. A more complete description of the STAR output files and the columns in them can be found in the STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

### Class input file
Class input file is a file created by the `run_class` step in SICILIAN and contains the information for all chimeric and spliced alignments extracted from the STAR BAM file. The class input file is saved in both parquet and tsv formats (`class_input.tsv` and `class_input.pq`; same for `class_input_secondary` which contains only secondary alignments). Each row in the class input file contains the alignment positions, read id, and alignment features needed for statistical modeling for a spliced/chimeric alignment. Class input file is the input file for SICILIAN statistical modeling step activated by the `run_GLM` flag.  

**Junction classification**
There are four categories for junctions in SICILIAN, depending on the relative positions of acceptor and donor sites and how far apart they are from each other.   
* `linear`: acceptor and donor are on the same chromosome and strand, closer than 1 MB to each other, and are based on the reference genome canonical ordering 
* `rev`: acceptor and donor are on the same chromosome and strand, closer than 1 MB to each other, and are ordered opposite of the reference genome canonical ordering (evidence for circRNAs)
* `sc`: local strandcrosses in which acceptor and donor are on the same chromosome but opposite strands.
* `fusion`: fusion junctions in which acceptor and donor are on different chromosomes, or on the same chromosome and strand but farther than 1MB from each other. 

### GLM_output.txt: 
This file is built by the last step in SICILIAN which performs the actual statistical modeling to assign a statistical score to each candidate junction  (activated by `run_GLM` flag) and contains all junctions that exist in `class_input.tsv` along with their statistical scores and junction-level alignment summaries.
The following columns exist in the `GLM_output.txt` file:

* `refName_newR1`: the junction id that contains gene name, positions, and strands for both acceptor and donor sides of the junction
* `fileTypeR1`: equals `Aligned` if the junction comes from `Aligned` alignments in the BAM file and equals `Chimeric` if it comes from `Chimeric` alignments in the BAM file
* `juncPosR1A`: position of the donor (5') side of the junction 
* `juncPosR1B`: position of the acceptor (3') side of the junction
* `chrR1A`: chromosome of the donor (5') side of the junction
* `chrR1B`: chromosome of the acceptor (5') side of the junction
* `read_strandR1A`:
* `read_strandR1B`:
* `gene_strandR1A`: The strand that the gene at `juncPosR1A` is on (`+` or `-`) based on the GTF file
* `gene_strandR1B`: The strand that the gene at `juncPosR1B` is on (`+` or `-`) based on the GTF file
* `geneR1A_uniq`: gene name for the donor (5') side of the junction
* `geneR1B_uniq`: gene name for the acceptor (3') side of the junction
* `geneR1A_ensembl`: the gene ensembl id for `geneR1A_uniq`
* `geneR1B_ensembl`: the gene ensembl id for `geneR1B_uniq`
* `numReads`: number of reads that have been aligned to the junction
* `geneR1A_expression`: the gene counts (HTseq counts) for `geneR1A_uniq` according to column V3 in `1ReadsPerGene.out.tab` for PE data (`2ReadsPerGene.out.tab` for SE data)
* `geneR1B_expression`: the gene counts (HTseq counts) for `geneR1B_uniq` according to column V3 in `1ReadsPerGene.out.tab` for PE data (`2ReadsPerGene.out.tab` for SE data)
* `median_overlap_R1`: median of junction overlaps across reads aligned to the junction
* `sd_overlap`: the standard deviation of junction overlap across reads aligned to the junction
* `threeprime_partner_number_R1`: number of distinct 3' splice sites across the class input file for the 5' side of the junction
* `fiveprime_partner_number_R1`: number of distinct 5' splice sites across the class input file for the 3' side of the junction
* `p_predicted_glmnet_constrained`: aggregated score for the junction  
* `p_predicted_glmnet_corrected_constrained`: aggregated score for the junction with correction for anomalous reads (calculated only for paired-end data)
* `junc_cdf_glmnet_constrained`: cdf of `p_predicted_glmnet_constrained` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet_corrected_constrained`: cdf of `p_predicted_glmnet_corrected_constrained` relative to its null distribution by randomly assigning reads to junctions
* `frac_genomic_reads`: fraction of aligned reads to the junction that have been also mapped to a genomic region
* `ave_min_junc_14mer`: the average of `min_junc_14mer` across the reads aligned to the junction
* `ave_max_junc_14mer`: the average of `max_junc_14mer` across the reads aligned to the junction
* `ave_AT_run_R1`: the average of `AT_run_R1` across the reads aligned to the junction, where `AT_run_R1` is the length of the longest run of A's (or T's) for a read alignment in R1
* `ave_GC_run_R1`: the average of `GC_run_R1` across the reads aligned to the junction, where `GC_run_R1` is the length of the longest run of G's (or C's) for a read alignment in R1
* `ave_max_run_R1`: the average of `max_run_R1` or max(`AT_run_R1`, `GC_run_R1`) across the reads aligned to the junction in R1
* `ave_entropy_R1`: the average of read sequence entropy calculated based on 5-mers. 
* `frac_anomaly`: the fraction of the aligned reads for the junction that are anamolous
* `p_val_median_overlap_R1`: the p-value of the statistical test for comparing the median overlaps of the aligned reads to the junction against the null of randomly aligned reads and small p-values are desired as they indicate that the median_overlap of the junction is large enough.  
* `emp.p_glmnet_constrained`: tme empirical p-value obtained based on `junc_cdf_glmnet_constrained` (the statistical score used for calling junctions in SE data)
* `emp.p_glmnet_corrected_constrained`: tme empirical p-value obtained obtained based on `junc_cdf_glmnet_corrected_constrained` (the statistical score used for calling junctions in SE data)
* `seqR1`: a representative sequence for the junction (the sequence of one of the reads aligned to the junction)

### sicilian_called_splice_juncs.tsv: 
`sicilian_called_splice_juncs.tsv` is the final output file by SICILIAN that contains the splice junctions called based upon SICILIAN statistical modeling. This file also includes extra annotation information (for splice junction, exon boundaries, and protein domain) for the called junctions. In addition to the columns that are also present in `GLM_output.txt`, `sicilian_called_splice_juncs.tsv` contains the following columns for the annotation of the called junctions:   

* `splice_ann`: both sides of the junction are annotated as exon boundaries and the junction itself is found in the GTF; subset of `both_ann`
* `both_ann`:  both sides of the junction are annotated as exon boundaries; this is equivalent to `exon_annR1A` AND `exon_annR1B`
* `exon_annR1A`: the exon on the first half of the junction is at an annotated boundary
* `exon_annR1B`: the exon on the second half of the junction is at an annotated boundary
* `missing_domains`: the protein domains that have been missed (spliced out) due to splicing   
* `domain_insertions`: the protein domains to which the splice junction adds amino acid sequence
