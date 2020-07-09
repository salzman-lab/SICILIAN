# SICILIAN


This project aligns fastq files to the genome using STAR, concatenates the `SJ.out.tab` and `Chimeric.out.junction` outputs into one file, and creates a class input file based on the STAR output.

## How to run the script

The script `write_jobs.py` is the main wrapper that runs all necessary jobs. To run on a sample, edit the following variables in the main function:

### Inputting sample data
Update the script with information about your sample by inputting the following parameters.
* `data_path`: set equal to the path containing the fastqs. Example: `data_path = "/scratch/PI/horence/JuliaO/single_cell/data/SRA/19.05.31.GSE109774/"`
* `assembly`: set equal to the keyword for the genome assembly you want to use (maybe "mm10" for mouse assembly 10, or hg38 for human assembly 38). Example: `assembly = "mm10"`
* `gtf_file`: The path to the gtf file that should be used to annotate genes. Example: `gtf_file = "/share/PI/horence/circularRNApipeline_Cluster/index/mm10_genes.gtf"`
* `run_name`: The unique name you want to give to the current run. It can be useful to include the date and a signifier for the data. Example: `run_name = "GSE109774_colon"`
* `r_ends`: list of unique endings for the file names of read one (location 0) and read 2 (location 1). Example: `r_ends = ["_1.fastq.gz", "_2.fastq.gz"]`
* `names`: list of unique identifiers for each fastq; e.g. file location for read 1 should be <data_path><name><r_ends[0]> and file location for read 2 should be <data_path><name><r_ends[1]>. Example: `names = ["SRR65462{}".format(i) for i in range(73,74)]`
* `single`: Set equal to True if the data is single-end, and False if it is paired-end. Note that currently if `single = True` it is assumed that the single read to be aligned is in the second fastq file (because of the tendancy for droplet protocols to have read 1 contain the barcode and UMI and read 2 contain the sequence that aligns to the genome). This also causes the files to be demultiplexed to create a new fastq file before they're mapped.
    
### Choosing STAR parameters
These parameters modify which combinations of STAR parameters are run. Be careful about including too many; even if you just have two values for each, you will run the pipeline 16 times. If you have 12 samples, then you're running the pipeline 192 times:
* `chimSegmentMin`: A list containing every value of the `--chimSegmentMin` STAR parameter you want to run on. Example: `chimSegmentMin = [12,10]`
* `chimJunctionOverhangMin`: A list containing every value of the `--chimJunctionOverhangMin` STAR parameter you want to run on. Example: `chimJunctionOverhangMin = [13, 10]`
* `alignSJstitchMismatchNmax`: For every value `<x>` in this list, STAR will be run with `--alignSJstitchMismatchNmax <x> -1 <x> <x>`. Example: `alignSJstitchMismatchNmax = [0,1]`
* `chimSegmentReadGapMax`: A list containing every value of the `--chimSegmentReadGapMax` STAR parameter you want to run on. Example: `chimSegmentReadGapMax = [0,3]`

### Choosing which scripts to run
These parameters let you decide which portion of the script you want to run (for example, if you're modifying the `create_class_input.py` script only, so the mapping shouldn't change):
* `run_whitelist`: Set equal to True if you want to run UMI-tools whitelist script to extract cell barcodes and identify the most likely true cell barcodes (will be run only for 10X)
* `run_extract`: Set equal to True if you want to run UMI-tools extract script which removes UMIs from fastq reads and append them to read name
* `run_map`: Set equal to True if you want to run the mapping job, and False otherwise
* `run_star_fusion`: Set equal to True if you want to run the STAR-Fusion, and False otherwise
* `run_ann`: Set equal to True if you want to annotate and concatenate the STAR files, False otherwise
* `run_class`: Set equal to True if you want to create the class input file, false otherwise
* `run_modify_class`: Set equal to True if you want to make consistent junction IDs in which reverse strands and discrepancies between reporting alignments in chimeric and aligned sam files are taken care of.
* `run_ensembl`: Set equal to True if you want to add gene names to the STAR gene count file and add gene ensembl ids and gene counts to the class input file, False otherwise
* `run_compare`: Set equal to True if you want to comapre the junction in the class inpout file with those in the STAR and STAR-Fusion ouput files, false otherwise
* `run_GLM`: Set equal to True if you want to run the GLM step and assign statistical scores to each junction in the class input file. The output of this step is a file named `GLM_output.txt`. 

After assigning these variables, run `python3 write_jobs.py` to submit the jobs.

## Description of output

There will be a separate folder for every combination of STAR parameters that was run. These folders will follow the naming convention `output/<run_name>_cSM_<chimSegmentMin value>_cJOM_<chimJunctionOverhangMin value>_aSJMN_<alignSJstitchMismatchNmax value>_cSRGM_<chimSegmentReadGapMax value>`. Example: `output/GSE109774_colon_cSM_10_cJOM_10_aSJMN_1_cSRGM_3`. Each of these folders will contain a folder for each variable in `names`, so a different folder for each pair of fastq files. Each of these sub-folders contains the following:

### STAR output
Based on the parameters we are running with, this includes  `2Aligned.out.sam`, `2Chimeric.out.sam`, `2Chimeric.out.junction`, `2SJ.out.tab`, and the same files with 1 instead of 2 if the reads are paired-end (among other STAR-generated files).

### Concatenated STAR splice files
The created files are called `1SJ_Chimeric.out` (if the data is paired-end) and `2SJ_Chimeric.out`. Each of these concatenates the respective `SJ.out.tab` and `Chimeric.out.junction` files, and adds columns for the gene names of the donor and acceptor as well as their strands (because we are defining strand to be whatever strand the gene is on, and `?` if there is no gene on either strand or a gene on both strands). However, these two files don't have the same columns, so right now the only columns that are shared are "donor_chromosome", "acceptor_chromosome", "donor_gene", and "acceptor_gene". The rest of the columns from the individual file are still present, but they are left blank for rows that belong to the "other" file. A more complete description of the columns can be found in the STAR manual: http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf

### Class input file
The purpose of the class input files is for each row to contain a summary of read 1 and read 2 for every read where read 1 is either Chimeric or contains a gap. 

### STAR-Fusion files
STAR-Fusion is an additional module developed by STAR authors to call fusion junctions using the Chimeric.out.junction file generated by STAR. All STAR-Fusion output files are written in the star-fusion folder which is placed in the same folder that other STAR output files are written. The final fusion calls by STAR-Fusion can be found in star-fusion.fusion_predictions.abridged.tsv. Currently STAR-Fusion is only run for R1 (even for PE data).

<!---#### Aligned vs Chimeric priority
Each sample folder will contain a filed called `class_input_priorityAlign.tsv` and `class_input_priorityChimeric.tsv`. The only difference between these two files is that in `class_input_priorityAlign.tsv` if a read appears in both `Aligned.out.sam` and `Chimeric.out.sam` the information will be taken from the `Aligned` file; for `class_input_priorityChimeric.tsv` the information will be taken from the `Chimeric` file before the `Aligned` file (this goes for read 1 and read 2). --->

#### Determining whether a read is included in the class input file
For single-end reads, a read is included in the class input file if it is Chimeric (contains the `ch` flag in the BAM file) or has an N in its cigar string (N codes for a base being skipped). For paired-end reads, if read1 is Chimeric or has an N in its cigar string, the read will be included in the class input file. Its line in the class input file will also include information about read 2, regardless of whether r2 is Chimeric/has an N in its CIGAR string or not. Note that if read 1 is not Chimeric and doesn't have an N in its CIGAR string, then read 2 won't be included in the class input file **even if read 2 is Chimeric or has an N in its CIGAR string** (this is something we could change).

Each read is included in `class_input.tsv` at most one time. If a read has a genomic read and a junctional read, the junctional read will be included in the class input file (and the fact that a genomic alignment exists will be marked by a zero in the `genomicAlignmentR1` column). However, if there are multiple junctional reads then the read with the lowest value in the `HI` tag will be included in `class_input.tsv`. All other junctional alignments can be found in `class_input_secondary.tsv`.

<!---For single-end reads, a read is included in the class input file if it is in `Chimeric.out.sam` or if it has an N in its CIGAR string in `Aligned.out.sam` (N codes for a base being skipped). Note that if the read is in both `Aligned.out.sam` and its CIGAR string doesn't have an N, and `Chimeric.out.sam`, then the read won't be included in `class_input_priorityAlign.tsv` but will be included in `class_input_priorityChimerc.tsv`. For paired-end reads, if read 1 is in `1Chimeric.out.sam` or if it has an N in its CIGAR string, the read will be included in the class input file (again, if it also appears in `1Aligned.out.sam` without an N then it won't appear in `class_input_priorityAlign.tsv`). Its line in the class input file will also include information about read 2, regardless of whether r2 is Chimeric/has an N in its CIGAR string or not. Note that if read 1 is not in `1Chimeric.out.sam` and doesn't have an N in its CIGAR string in `1Aligned.out.sam`, the read won't be included in the class input file **even if read 2 is Chimeric or has an N in its CIGAR string** (this is something we could change). ---> 

#### Fields of the class input file

The class input file is saved in both parquet format and tsv format (`class_input.tsv` and `class_input.pq`; same for `class_input_secondary`). `class_input.tsv` is modified by the GLM script to include the SICILIAN score, to deduplicate UMIs (by UMI + barcode + junction name), and to add a few columns such as `overlap_R1`. However, `class_input.pq` is not modified by the GLM.

A note on the naming convention for the fields of the class input file: `id` and `class` are the only fields that are necessarily the same for read 1 and read 2. All other fields have `R1` in them if they pertain to read 1, and `R2` in them if they pertain to read 2. For single-end reads, the information for the one read will always show up in the `R1` columns, even if it's actually from the fastq file labeled 2. Then within read 1 and read 2 the columns are split into `A` and `B`. For a read that aligns to two locations (either in the Chimeric file, or with an N in the CIGAR string in the Aligned file), the first portion of the read is referred to as `A`  and the second is referred to as `B`. Here "first portion" means that if you saw the read in the raw fastq file, the first bases in the read would align to the `A` location, and the last bases would align to the `B` location.

Categories of the form `<field>R2` rather than `<field>R1` follow the same definition but for the second read unless otherwise indicated.

If a read is from the Aligned file rather than the Chimeric file, the following columns will have the value `NA`: `flagB`, `posB`, `aScoreB`, `qualB`, `seqB`. If a read doesn't contain a junction, the following columns will **also** have the value `NA`: `strandR2B`, `chrR2B`, `geneR2B`, `readClassR2`, `juncPosR2A`, `juncPosR2B`, `cigarR2B`, `MR2B`, `SR2B`. 

The class input file contains the following fields (not necessarily in this order):

<!---Note: :whale: indicates that it is safe to use without modification in a model between aligned and chimeric. :octopus: means it is safe to use without modification within chimeric. ---> 

**Junction classification**
* `linear`: acceptor and donor are on the same chromosome and strand, closer than 1 MB to each other, and are based on the reference genome canonical ordering 
* `rev`: acceptor and donor are on the same chromosome and strand, closer than 1 MB to each other, and are ordered opposite of the reference genome canonical ordering
* `sc`: acceptor and donor are on the same chromosome but opposite strands.
* `fusion`: acceptor and donor are on different chromosomes, or on the same chromosome and strand but farther than 1MB from each other. 

**Refnames for negative-strand genes will have the acceptor first and the donor second**
* `id`: Read name. Example: `SRR6546273.367739`
* `readLenR1`: length of read 1 (including any softclipped portions)
* `fileTypeR1`: equals `Aligned` if the read came from the `1Aligned.out.sam` file, `Chimeric` if it came from the `1Chimeric.out.sam` file.
* `seqR1`: The read sequence
* `AT_run_R1`: max(the length of the longest run of A's, the length of the longest run of T's). The A's and T's are combined because we could have the sequence or the reverse complement. This is calculated from `seqR1`.
* `GC_run_R1`: max(the length of the longest run of G's, the length of the longest run of C's). This is calculated from `seqR1`
* `max_run_R1`: max(`AT_run_R1`, `GC_run_R1`)
* `entropyR1`: The entropy of the read calculated based on 5-mers. Let $k_1, \ldots, k_n$ be all the 5-mers in the read sequence (for example, for ACTCCGAGTCCTCCG the 5-mers would be **ACTCC**GAGTCCTCCG, A**CTCCG**AGTCCTCCG, AC**TCCGA**GTCCTCCG, ACT**CCGAG**TCCTCCG, ACTC**CGAGT**CCTCCG, ACTCC**GAGTC**CTCCG, ACTCCG**AGTCC**TCCG, ACTCCGA**GTCCT**CCG, ACTCCGAG**TCCTC**CG, ACTCCGAGT**CCTCC**G, and ACTCCGAGTC**CTCCG**). Then let $N(c)$ equal the number of times the kmer c appears in $k_1, \ldots, k_n$ (for example, CTCCG appears twice in ACTCCGAGTCCTCCG). Then the entropy is defined as $\sum_C -\frac{N(c)}{n}\log\left(\frac{N(c)}{n}\right)$ (here C is the set of unique 5-mers in the read)
* `refName_ABR1`: The refName for R1 will always be of the form `<chrR1A>:<geneR1A>:<juncPosR1A>:<gene_strandR1A>|<chrR1B>:<geneR1B>:<juncPosR1B>:<gene_strandR1B>|<readClassR1>`
* `UMI`: for 10X data this is the UMI found through UMI-tools, otherwise this is NA
* `barcode`: for 10X data this is the barcode found through UMI-tools, otherwise this is NA
* `maxA_10merR1`: The maximum number of A's in a 10-base stretch in the read
* `maxT_10merR1`: The maximum number of T's in a 10-base stretch in the read
* `maxG_10merR1`: The maximum number of G's in a 10-base stretch in the read
* `maxC_10merR1`: The maximum number of C's in a 10-base stretch in the read
* `aScoreR1A`: alignment score from the SAM file after the `AS`.
* `aScoreR1B`: alignment score from the SAM file after the `AS`.
* `MR1A`: The number of M's in `cigarR1A`
* `MR1B`: The number of M's in `cigarR1B`
* `SR1A`: The number of S's in the portion of the cigar string corresponding to the first half of the splice
* `SR1B`: The number of S's in the portion of the cigar string corresponding to the second half of the splice
* `nmmR1A`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
* `nmmR1B`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
* `qualR1A`: Mapping quality of the first portion of read 1. From the manual: "The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10\*log10(1-1/Nmap)) for multi-mapping reads" 
* `qualR1B`: Mapping quality for the second portion of read 1.
* `NHR1A`: Number of reported alignments that contains the query in the current record
* `NHR1B`: Number of reported alignments that contains the query in the current record
* `HIR1A`: Query hit index, indicating the alignment record is the i-th one stored in SAM
* `HIR1B`: Query hit index, indicating the alignment record is the i-th one stored in SAM
* `cigarR1A`: The cigar string for portion A (for Chimeric, this is without the softclipped portion that corresponds to B; for Aligned, this is without the long N sequence marking the intron and everything after)
* `cigarR2B`: The cigar string for portion B
* `juncPosR1A`: The last position part A of read 1 aligns to before the junction. If `fileTypeR1 = Chimeric`: if `flagR1A` is 0 or 256, this is equal to `posR1A + ` the sum of the M's, N's, and D's in the CIGAR string. If `flagR1A` is 16 or 272 this is equal to `posR1A`. If `fileTypeR1 = Aligned`, then this equals `posR1A` plus the sum of the M's, N's, and D's before the largest N value in the CIGAR string. 
* `juncPosR1B`: The first position of part B of read 1 that aligns after the junction. If `fileTypeR1 = Chimeric`: if `flagR1B` is 0 or 256, this is equal to `posR1B`. If `flagR1B` is 16 or 272 this is equal to `posR1B + ` the sum of the M's, N's, and D's in the CIGAR string. If `fileTypeR1 = Aligned`, then this equals `posR1B`. 
* `geneR1A`: Gene that read 1 part A was aligned to. If no gene was annotated in that area, it's marked as "unknown". If multiple genes are annotated in this area, it's marked with all of those gene names concatenated with commas in between Example: `Ubb,Gm1821`. Also see the note on the annotation. 
* `geneR1B`: Gene that read 1 part B was aligned to.
* `chrR1A`: Chromosome that read 1 part A was aligned to
* `chrR1B`: Chromosome that read 1 part B was aligned to
* `read_strand_compatible`: A 1 here indicates that the read strands are compatible, 0 indicates that they're not. This is 1 if `read_strandR1A` doesn't equal `read_strandR2A` and 0 otherwise. Note that this is **only** based on the "A" part of the read; cross-reference with `strand_crossR1` and `strand_crossR2` to verify that all read strands are in accordance. This is not present for single-end reads.
* `location_compatible`: A 1 here indicates that the locations are compatible, a 0 indicates that they're not. If `readClassR1` is `fus`, this is 1. If `readClassR1` is `sc`, this is 0. If `readClassR1` is `lin` and `read_strandR1` is +, then this is 1 if `posR1A` < `posR2A` and 0 otherwise. If `readClassR1` is `lin` and `read_strandR1` is -, then this is 1 if `posR1A` > `posR2A` and 0 otherwise. if `readClassR1` is `rev`, then this is 1 if `posR2A` is between `posR1A` and `posR1B` and 0 otherwise. Note that this is **only** based on the "A" part of read 2 for right now. This is not present for single-end reads.
* `genomicAlignmentR1`: Here 1 indicates that there is a genomic alignment for read 1, and 0 indicates that there is not.
* `primaryR1A`: 
* `primaryR1B`:
* `genomic_aScoreR1`: the maximum alignment score of all genomic alignments of that read (NA if `genomicAlignmentR1` == 0)
* `spliceDist`: The absolute difference between `juncPosR1A` and `juncPosR1B`
* `refName_newR1`: the consistent disambiguated junction id obtained by `run_modify_class` step. All subsequent analyses on junctions are based on the ids in this column. 
* `gene_strandR1A`: The strand that the gene at `juncPosR1A` is on (`+` or `-`); if there is no gene at that location, or there is a gene on both strands, this equals `?`.
* `gene_strandR1B`: The strand that the gene at `juncPosR1B` is on; if there is no gene at that location, or there is a gene on both strands, this equals `?`.
* `geneR1A_uniq`: Gene name after disambiguation steps have been run on it (to try to narrow down to just one gene)
* `geneR1B_uniq`: Gene name after disambiguation steps have been run on it (to try to narrow down to just one gene)
* `max_id_priority`: The minimum `HIR1A` value for a given `id` (used to split alignments between `class_input` and `class_input_secondary`)
<!---
1. `geneR1B_uniq`: 
2. `geneR1A_uniq`:
3. 
4. `class`: **deprecated: see read_strand_compatible and location_compatible** Class defined by read 1 or read 2. For paired end mode, the options are circular, linear, decoy, err (this happens when the strand is ambiguous, because in that case we can't tell if a potential circular junction is a circle or a decoy, since this definition depends on the strand), fusion (read 1 and read 2 are on different chromosomes, or either r1 or r2 is split between two chromosomes), and strandCross (both reads have flags indicating they're both on + or both on -; this is before we correct strandedness by gene location). For single end data the options are lin (linear-type junction), rev (circle-type junction), and fus (part of read maps to one chromosome and part maps to another)
5. `refName_ABR1`: The refName for R1 will always be of the form `<chrR1A>:<geneR1A>:<juncPosR1A>:<gene_strandR1A>|<chrR1B>:<geneR1B>:<juncPosR1B>:<gene_strandR1B>|<readClassR1>`
6. `refName_readStrandR1`: The refName for R1 will always be of the form either `<chrR1A>:<geneR1A>:<juncPosR1A>:<gene_strandR1A>|<chrR1B>:<geneR1B>:<juncPosR1B>:<gene_strandR1B>|<readClassR1>` or `<chrR1B>:<geneR1B>:<juncPosR1B>:<gene_strandR1B>|<chrR1A>:<geneR1A>:<juncPosR1A>:<gene_strandR1A>|<readClassR1>`. Which one of these two it is will be defined by whether the read strand is + (in which case it will be the first one) or the read strand is - (in which case it will be the second one). See the descriptions for these individual columns for more specifics. **Refnames for negative-strand genes will have the acceptor first and the donor second**
7. `refName_ABR2`: The refName for R2 will always be of the form  `<chrR2A>:<geneR2A>:<juncPosR2A>:<gene_strandR2A>|<chrR2B>:<geneR2B>:<juncPosR2B>:<gene_strandR2B>|<readClassR2>`
8. `refName_readStrandR2`: If read 2 is from `2Chimeric.out.sam` or has an N in the CIGAR string, it will have the same format as `refNameR1` except `R2` replaces `R1`. If instead it aligns without gaps/chimera, it's name is `<chrR2A>:<geneR2A>:<gene_strandR2A>`.
9. `fileTypeR1`: equals `Aligned` if the read came from the `1Aligned.out.sam` file, `Chimeric` if it came from the `1Chimeric.out.sam` file.
10. `fileTypeR2`: equals `Aligned` if the read came from the `2Aligned.out.sam` file, `Chimeric` if it came from the `2Chimeric.out.sam` file.
11 `readClassR1`: This is `fus` if read 1 part A and read 1 part B map to different chromosomes. It is `sc`  if the flags of read 1 part A and read 1 part B don't match (they match if they're both in [0,256] or they're both in [16,276]; they don't match otherwise). It is `rev` if the positions are consistent with a circular junction. It is `lin` if the positions are consistent with a linear junction. It is `err` otherwise (this can occur due to `?` occuring as the strand)
12. `readClassR2`: See `readClassR1` and replace read 1 with read 2. If read 2 isn't a junction this will be `NA`.
13. `numNR1`: Number of N in the reference; from the XN tag in the SAM file
14. `numNR2`: Number of N in the reference; from the XN tag in the SAM file
15. `readLenR1`: length of read 1 (including any softclipped portions)
16. `readLenR2`: Length of read 2 (including any softclipped portions)
17. `barcode`: for 10X data this is the barcode found through UMI-tools, otherwise this is NA
18. `UMI`: for 10X data this is the UMI found through UMI-tools, otherwise this is NA
19. `entropyR1`: The entropy of the read calculated based on 5-mers. Let $k_1, \ldots, k_n$ be all the 5-mers in the read sequence (for example, for ACTCCGAGTCCTCCG the 5-mers would be **ACTCC**GAGTCCTCCG, A**CTCCG**AGTCCTCCG, AC**TCCGA**GTCCTCCG, ACT**CCGAG**TCCTCCG, ACTC**CGAGT**CCTCCG, ACTCC**GAGTC**CTCCG, ACTCCG**AGTCC**TCCG, ACTCCGA**GTCCT**CCG, ACTCCGAG**TCCTC**CG, ACTCCGAGT**CCTCC**G, and ACTCCGAGTC**CTCCG**). Then let $N(c)$ equal the number of times the kmer c appears in $k_1, \ldots, k_n$ (for example, CTCCG appears twice in ACTCCGAGTCCTCCG). Then the entropy is defined as $\sum_C -\frac{N(c)}{n}\log\left(\frac{N(c)}{n}\right)$ (here C is the set of unique 5-mers in the read)
20. `entropyR2`: Same as above
21. `seqR1`: The read sequence
22. `seqR2`: The read sequence
23. `read_strand_compatible`: A 1 here indicates that the read strands are compatible, 0 indicates that they're not. This is 1 if `read_strandR1A` doesn't equal `read_strandR2A` and 0 otherwise. Note that this is **only** based on the "A" part of the read; cross-reference with `strand_crossR1` and `strand_crossR2` to verify that all read strands are in accordance. This is NA for single-end reads.
24. `location_compatible`: A 1 here indicates that the locations are compatible, a 0 indicates that they're not. If `readClassR1` is `fus`, this is 1. If `readClassR1` is `sc`, this is 0. If `readClassR1` is `lin` and `read_strandR1` is +, then this is 1 if `posR1A` < `posR2A` and 0 otherwise. If `readClassR1` is `lin` and `read_strandR1` is -, then this is 1 if `posR1A` > `posR2A` and 0 otherwise. if `readClassR1` is `rev`, then this is 1 if `posR2A` is between `posR1A` and `posR1B` and 0 otherwise. Note that this is **only** based on the "A" part of read 2 for right now. This is NA for single-end reads.
25. `strand_crossR1`: Here 1 indicates that there was a strand cross in read 1 (meaning `read_strandR1A` doesn't equal `read_strandR1B`) and 0 indicates that there wasn't a strand cross (meaning the read strands are equal) 
26. `strand_crossR2`: Here 1 indicates that there was a strand cross in read 1 (meaning `read_strandR2A` doesn't equal `read_strandR2B`) and 0 indicates that there wasn't a strand cross (meaning the read strands are equal)
27. `genomicAlignmentR1`: Here 1 indicates that there is a genomic alignment for read 1, and 0 indicates that there is not. 
28. `spliceDist`: the absolute value of the difference between `juncPosR1A` and `juncPosR1B` 
29. `AT_run_R1`: max(the length of the longest run of A's, the length of the longest run of T's). The A's and T's are combined because we could have the sequence or the reverse complement. This is calculated from `seqR1`.
30. `GC_run_R1`: max(the length of the longest run of G's, the length of the longest run of C's). This is calculated from `seqR1`
31. `max_run_R1`: max(`AT_run_R1`, `GC_run_R1`)
32. `AT_run_R2`:  max(the length of the longest run of A's, the length of the longest run of T's). The A's and T's are combined because we could have the sequence or the reverse complement. This is calculated from `seqR2`.
33. `GC_run_R2`: max(the length of the longest run of G's, the length of the longest run of C's). This is calculated from `seqR2`
34. `max_run_R2`: max(`AT_run_R2`, `GC_run_R2`)
35. `Organ`:
36. `Cell_Type.s.`:
37. `chrR1A`: Chromosome that read 1 part A was aligned to
38. `chrR1B`: Chromosome that read 1 part B was aligned to
39. `chrR2A`: Chromosome that read 2 part A was aligned to
40. `chrR2B`: Chromosome that read 2 part B was aligned to
41. `geneR1A`: Gene that read 1 part A was aligned to. If no gene was annotated in that area, it's marked as "unknown". If multiple genes are annotated in this area, it's marked with all of those gene names concatenated with commas in between Example: `Ubb,Gm1821`. Also see the note on the annotation. 
42. `geneR1B`: Gene that read 1 part B was aligned to.
43. `geneR2A`: Gene that read 2 part A was aligned to.
44. `geneR2B`: Gene that read 2 part B was aligned to.
45. `juncPosR1A`: The last position part A of read 1 aligns to before the junction. If `fileTypeR1 = Chimeric`: if `flagR1A` is 0 or 256, this is equal to `posR1A + ` the sum of the M's, N's, and D's in the CIGAR string. If `flagR1A` is 16 or 272 this is equal to `posR1A`. If `fileTypeR1 = Aligned`, then this equals `posR1A` plus the sum of the M's, N's, and D's before the largest N value in the CIGAR string. 
46. `juncPosR1B`: The first position of part B of read 1 that aligns after the junction. If `fileTypeR1 = Chimeric`: if `flagR1B` is 0 or 256, this is equal to `posR1B`. If `flagR1B` is 16 or 272 this is equal to `posR1B + ` the sum of the M's, N's, and D's in the CIGAR string. If `fileTypeR1 = Aligned`, then this equals `posR1B`. 
47. `juncPosR2A`: This follows the same rules as those for read 1, except if the read doesn't contain a junction this is `NA`
48. `juncPosR2B`: This follows the same rules as those for read 1, except if the read doesn't contain a junction this is `NA`
49. `gene_strandR1A`: The strand that the gene at `juncPosR1A` is on (`+` or `-`); if there is no gene at that location, or there is a gene on both strands, this equals `?`.
50. `gene_strandR1B`: The strand that the gene at `juncPosR1B` is on; if there is no gene at that location, or there is a gene on both strands, this equals `?`.
51. `gene_strandR2A`: The strand that the gene at `juncPosR2A` is on; if there is no gene at that location, or there is a gene on both strands, this equals `?`.
52. `gene_strandR2B`: The strand that the gene at `juncPosR2B` is on; if there is no gene at that location, or there is a gene on both strands, this equals `?`.
53. :octopus: `aScoreR1A`: alignment score from the SAM file after the `AS`.
54. :octopus: `aScoreR1B`: alignment score from the SAM file after the `AS`.
55. `aScoreR2A`: alignment score from the SAM file after the `AS`.
56. `aScoreR2B`: alignment score from the SAM file after the `AS`.
57. `flagR1A`: flag from the SAM file; 0 means forward strand primary alignment, 256 means forward strand secondary alignment, 16 means reverse strand primary alignment, 272 means reverse strand secondary alignment.
58. `flagR1B`: flag from the SAM file; 0 means forward strand primary alignment, 256 means forward strand secondary alignment, 16 means reverse strand primary alignment, 272 means reverse strand secondary alignment.
59. `flagR2A`: flag from the SAM file; 0 means forward strand primary alignment, 256 means forward strand secondary alignment, 16 means reverse strand primary alignment, 272 means reverse strand secondary alignment.
60. `flagR2B`: flag from the SAM file; 0 means forward strand primary alignment, 256 means forward strand secondary alignment, 16 means reverse strand primary alignment, 272 means reverse strand secondary alignment.
61. `posR1A`: 1-based leftmost mapping position from the SAM file 
62. `posR1B`: 1-based leftmost mapping position from the SAM file 
63. `posR2A`: 1-based leftmost mapping position from the SAM file 
64. `posR2B`: 1-based leftmost mapping position from the SAM file 
65. `qualR1A`: Mapping quality of the first portion of read 1. From the manual: "The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10\*log10(1-1/Nmap)) for multi-mapping reads" 
66. `qualR1B`: Mapping quality for the second portion of read 1.
67. `qualR2A`: Mapping quality of the first portion of read 2. 
68. `qualR2B`: Mapping quality of the second portion of read 2.
69.`MDR1A`: The MD flag from the SAM file (indicates where mutations, insertions, and delections occur)
70. `MDR1B`: The MD flag from the SAM file (indicates where mutations, insertions, and delections occur)
71. `MDR2A`: The MD flag from the SAM file (indicates where mutations, insertions, and delections occur)
72. `MDR2B`: The MD flag from the SAM file (indicates where mutations, insertions, and delections occur)
73. :octopus: `nmmR1A`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
74. :octopus: `nmmR1B`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
75. `nmmR2A`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
76. `nmmR2B`: The number of mismatches in the read; calculated by finding the number of times A,C,T or G appears in the MD flag
77. `cigarR1A`: The cigar string for portion A (for Chimeric, this is without the softclipped portion that corresponds to B; for Aligned, this is without the long N sequence marking the intron and everything after)
78. `cigarR2B`: The cigar string for portion B
79. `cigarR2A`: The cigar string for portion A
80. `cigarR2B`: The cigar string for portion B
81. :whale: :octopus: `MR1A`: The number of M's in `cigarR1A` (this corresponds to the number of bases that have a match or mismatch with the reference)
82. :whale: :octopus: `MR1B`: The number of M's in `cigarR1B`
83. `MR2A`: The number of M's in `cigarR2A`
84. `MR2B`: The number of M's in `cigarR2B`
85. :whale: :octopus: `SR1A`: The number of S's in `cigarR1A` (this corresponds to the number of bases that have been softclipped)
86. :whale: :octopus: `SR1B`: The number of S's in `cigarR1B`
87. `SR2A`: The number of S's in `cigarR2A`
88. `SR2B`: The number of S's in `cigarR2B`
89. `NHR1A`: Number of reported alignments that contains the query in the current record
90. `NHR1B`: Number of reported alignments that contains the query in the current record
91. `NHR2A`: Number of reported alignments that contains the query in the current record
92. `NHR2B`: Number of reported alignments that contains the query in the current record
93. `HIR1A`: Query hit index, indicating the alignment record is the i-th one stored in SAM
94. `HIR1B`: Query hit index, indicating the alignment record is the i-th one stored in SAM
95. `HIR2A`: Query hit index, indicating the alignment record is the i-th one stored in SAM
96. `HIR2B`: Query hit index, indicating the alignment record is the i-th one stored in SAM
97. :octopus: `nMR1A`: The number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate
98. :octopus: `nMR1B`: The number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate
99. `nMR2A`: The number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate
100. `nMR2B`: The number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate
101. :octopus: `NMR1A`: Edit distance to the reference, including ambiguous bases but excluding clipping
102. :octopus: `NMR1B`: Edit distance to the reference, including ambiguous bases but excluding clipping
103. `NMR2A`: Edit distance to the reference, including ambiguous bases but excluding clipping
104. `NMR2B`: Edit distance to the reference, including ambiguous bases but excluding clipping
105. :octopus: `jMR1A`: intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
jI:B:I,Start2,End1,Start2,End2,...
106. :octopus: `jMR1B`: intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
jI:B:I,Start1,End1,Start2,End2,...
107. `jMR2A`: intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
jI:B:I,Start1,End1,Start2,End2,...
108. `jMR2B`: intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
jI:B:I,Start1,End1,Start2,End2,...
109. :octopus: `jIR1A`: attributes require samtools 0.1.18 or later, and were reported to be incompatible with some downstream tools such as Cufflink
110. :octopus: `jIR1B`: attributes require samtools 0.1.18 or later, and were reported to be incompatible with some downstream tools such as Cufflink
111. `jIR2A`: attributes require samtools 0.1.18 or later, and were reported to be incompatible with some downstream tools such as Cufflink
112. `jIR2B`: attributes require samtools 0.1.18 or later, and were reported to be incompatible with some downstream tools such as Cufflink
113. `read_strandR1A`: The read strand
114. `read_strandR1B`: The read strand
115. `read_strandR2A`: The read strand
116. `read_strandR2B`: The read strand
117. `gene_strandR1A_new`: the disambiguated strand for gene R1A (obtained by `run_modify_class` step) 
118. `gene_strandR1B_new`: the disambiguated strand for gene R1B (obtained by `run_modify_class` step)
119. `refName_newR1`: the consistent disambiguated junction id obtained by `run_modify_class` step. All subsequent analyses on junctions are based on the ids in this column. 
120. `is.STAR_Chim`: is 1 if a chimeric junction in the class input file is also present in `Chimeric.out.junction`, and is 0 otherwise. (It is NA for all splice alignments in the class inpout file coming from `Aligned.out.sam`)
121. `is.STAR_SJ`: is 1 if a splice junction in the class input file is also present in `SJ.out.tab`, and is 0 otherwise. (It is NA for all * chimeric junctions in the class inpout file coming from `Chimeric.out.sam`)
122. `intron_motif`: intron motif (obtained from `SJ.out.tab`) 
123. `is.annotated`: indicates whether the junction is annotated (obtained from `SJ.out.tab`) 
124. `num_uniq_map_reads`: number of uniquely mapping reads (obtained from `SJ.out.tab`)
125. `num_multi_map_reads`: number of multimapping reads (obtained from `SJ.out.tab`)
126. `maximum_SJ_overhang`: maximum overlap among aligned reads (obtained from `SJ.out.tab`)
127. `is.STAR_Fusion`: is 1 if a chimeric junction in the class input file is also present in `star-fusion.fusion_predictions.abridged.tsv` (called by STAR-Fusion), and is 0 otherwise. (It is NA for all splice alignments in the class inpout file coming from `Aligned.out.sam`)
128. `numReads`: number of reads aligned to the junction based on collapsing junction ids in 
129. `length_adj_AS_R1`: length adjusted alignment score for R1 
130. `length_adj_AS_R1A`: length adjusted alignment score for R1A
131. `length_adj_AS_R1B`: length adjusted alignment score for R1B
132. `nmmR1`: number of mismatches in R1
133. `qual_R1`: mapping quality value for R1
--->
* `overlap_R1`: junction overlap for R1 (the minimum of junction anchors)
* `max_overlap_R1`: maximum overlap for R1 (the maximum of junction anchors)
* `median_overlap_R1`: median of `overlap_R1` across all aligned reads for each junction 
* `is.zero_nmm`: indicates whether the R1 alignment is mismatch free
* `is.multimapping`: indicates whether the alignment is multimapping 
* `njunc_binR1A`: junction noise score for `R1A`
* `njunc_binR1B`: junction noise score for `R1B`
* `threeprime_partner_number_R1`: number of distinct 3' splice sites across the class input file for each 5' splice site in the class input file
* `fiveprime_partner_number_R1`: number of distinct 5' splice sites across the class input file for each 3' splice site in the class input file
* `length_adj_AS_R2`: length ajusted alignment score for R2
* `glm_per_read_prob`: per-read score for each read alignment computed by GLM 
* `glm_per_read_prob_corrected`: per-read score for each read alignment computed by GLM where per-read scores for anomalous reads are downscaled   
* `glmnet_per_read_prob`: per-read score for each read alignment computed by GLMnet 
* `glmnet_per_read_prob_corrected`: per-read score for each read alignment computed by GLMnet and per-read scores for anomalous reads are downscaled
* `glmnet_per_read_prob_constrained`: per-read score for each read alignment computed by constrained GLMnet model
* `glmnet_per_read_prob_corrected_constrained`: per-read score for each read alignment computed by constrained GLMnet model and per-read scores for anomalous reads are downscaled
* `glmnet_twostep_per_read_prob`: the two-step per-read score for each chimeric read 
* `glmnet_twostep_per_read_prob_constrained`:  per-read score for each read alignment computed by twostep GLMnet with constrained cofficients
* `p_predicted_glm`: aggregated score for each junction based on `glm_per_read_prob` 
* `p_predicted_glm_corrected`: aggregated score for each junction based on `glm_per_read_prob_corrected`
* `p_predicted_glmnet`: aggregated score for each junction based on `glmnet_per_read_prob`
* `p_predicted_glmnet_corrected`: aggregated score for each junction based on `glmnet_per_read_prob_corrected`
* `p_predicted_glmnet_constrained`: aggregated score for each junction based on `glmnet_per_read_prob_constrained`
* `p_predicted_glmnet_corrected_constrained`: aggregated score for each junction based on `glmnet_per_read_prob_corrected_constrained`
* `p_predicted_glmnet_twostep`: aggregated score for each junction based on `glmnet_twostep_per_read_prob`
* `p_predicted_glmnet_twostep_constrained`:  aggregated score for each junction based on `glmnet_twostep_per_read_prob_constrained`
* `junc_cdf_glm`: cdf of `p_predicted_glm` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glm_corrected`: cdf of `p_predicted_glm_corrected` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet`: cdf of `p_predicted_glmnet` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet_corrected`: cdf of `p_predicted_glmnet_corrected` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet_constrained`: cdf of `p_predicted_glmnet_constrained` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet_corrected_constrained`: cdf of `p_predicted_glmnet_corrected_constrained` relative to its null distribution by randomly assigning reads to junctions
* `junc_cdf_glmnet_twostep`: cdf of `p_predicted_glmnet_twostep` relative to its null distribution by randomly assigning reads to junctions
* `genomic_aScoreR1`: the maximum alignment score of all genomic alignments of that read (NA if `genomicAlignmentR1` == 0)
* `frac_genomic_reads`: fraction of aligned reads for each junction that have also genomic alignment 
* `frac_anomaly`: the fraction of the aligned reads for the junction that are anamolous
* `ave_min_junc_14mer`:  the average of `min_junc_14mer` across the reads aligned to the junction
* `ave_max_junc_14mer`:  the average of `max_junc_14mer` across the reads aligned to the junction
* `ave_AT_run_R1`:  the average of `AT_run_R1` across the reads aligned to the junction
* `ave_GC_run_R1`:  the average of `GC_run_R1` across the reads aligned to the junction
* `ave_max_run_R1`: the average of `max_run_R1` across the reads aligned to the junction
* `p_val_median_overlap_R1`: the p-value of the statistical test for comparing the median overlaps of the aligned reads to the junction against the null of randomly aligned reads and small p-values are desired as they indicate that the median_overlap of the junction is large enough.  
* `uniformity_test_pval`: the p_value of the uniformity test for the junction overlap of the reads aligned to the junction computed by chisq.test. It is computed only for junctions with at most 15 reads and p-values close to 1 are desired as they indicate reads are uniformly distributed.
* `sd_overlap`: the standard deviation of the junction overlap for the reads aligned to the junction 


### New columns in the class input file after run_ensembl step:
* `geneR1B_ensembl`: the gene ensembl id for `geneR1B`
* `geneR1A_ensembl`: the gene ensembl id for `geneR1A`
* `geneR1B_name`: the HUGO name for `geneR1B`
* `geneR1A_name`: the HUGO name for `geneR1A`
* `geneR1A_expression`: the gene counts (htseq counts) for `geneR1A` according to column V3 in `1ReadsPerGene.out.tab` (`2ReadsPerGene.out.tab` for PE data)
* `geneR1B_expression`: the gene counts (htseq counts) for `geneR1B` according to column V3 in `1ReadsPerGene.out.tab` (`2ReadsPerGene.out.tab` for PE data)


### GLM_output.txt: 
This file is built by the `GLM_script.R` and contains all junction level metrics including aggregated p_predicted and junc_cdf scores and other alignment qualities at the junction level. Specifically `GLM_output.txt` contains emperical p-values computed by comparing the `junc_cdf` scores obtained via various varioations of the GLM(net) model against the distribution of those scores in "Bad junctios", junctions with at least 10% of reads with genomic alignment. There are currently the following empirical p-values in the GLM output file:
* `emp.p_glm`: obtained based on `junc_cdf_glm`
* `emp.p_glmnet`: obtained based on `junc_cdf_glmnet`
* `emp.p_glmnet_constrained`: obtained based on `junc_cdf_glmnet_constrained`
* `emp.p_glm_corrected`: obtained based on `junc_cdf_glm_corrected`
* `emp.p_glmnet_corrected`: obtained based on `junc_cdf_glmnet_corrected`
* `emp.p_glmnet_corrected_constrained`: obtained based on `junc_cdf_glmnet_corrected_constrained`

 The following columns from the `class input file` is written in `GLM_output.txt`:

`refName_newR1`, `geneR1B_uniq`, `geneR1A_uniq`, `geneR1B_ensembl`, `geneR1A_ensembl`, `readClassR1`, `geneR1A`, `geneR1B`, `geneR1A_expression_unstranded`, `geneR1A_expression_stranded`, `geneR1B_expression_unstranded`, `geneR1B_expression_stranded`, `is.STAR_Chim`, `is.STAR_SJ`, `is.STAR_Fusion`, `numReads`, `median_overlap_R1`, `sd_overlap`, `njunc_binR1A`, `njunc_binR1B`, `threeprime_partner_number_R1`, `fiveprime_partner_number_R1`,`p_predicted_glm`, `junc_cdf_glm`, `p_predicted_glm_corrected`, `junc_cdf_glm_corrected`, `p_predicted_glmnet`, `junc_cdf_glmnet`, ` p_predicted_glmnet_corrected`, `junc_cdf_glmnet_corrected`, `p_predicted_glmnet_twostep`, `junc_cdf_glmnet_twostep`, `frac_genomic_reads`, `frac_anomaly`, `ave_min_junc_14mer`, `ave_max_junc_14mer`, `ave_AT_run_R1`, `ave_GC_run_R1`, `ave_max_run_R1`, `ave_AT_run_R2`, `ave_GC_run_R2`, `ave_max_run_R2`, `p_val_median_overlap_R1`, `uniformity_test_pval`,  `p_predicted_glmnet_constrained`, `junc_cdf_glmnet_constrained`, `p_predicted_glmnet_corrected_constrained`, `junc_cdf_glmnet_corrected_constrained`, `p_predicted_glmnet_twostep_constrained`, `seqR1`, `seqR2`

The following columsn are added from the STAR output file for "Aligned" junctions (`SJ.out.tab`): 

`intron_motif`, `is.annotated`, `num_uniq_map_reads`, `num_multi_map_reads`, `maximum_SJ_overhang`,


###  Files for comparing the junctions in the class input files with those in the STAR output files:
If run_compare is set to TRUE, junctions in both class input files are comapred with the junctions in the STAR output files `SJ.out.tab` and `Chimeric.out.junction`. Currently, the comaprison is only for junctional R1 reads. This step adds 3 columns to the class input file: `is.STAR_Chim`, `is.STAR_SJ`, and `is.STAR_Fusion`.

Also, the following 4 files will be written at the end of this module:
* `in_star_chim_not_in_classinput_priority_Align.txt`: chimeric junctions in the STAR `Chimeric.out.junction` file that cannot be found in `class_input_priorityAlign.tsv`
* `in_star_chim_not_in_classinput_priority_Chim.txt`: chimeric junctions in the STAR `Chimeric.out.junction` file that cannot be found in `class_input_priorityChimeric.tsv`
* `in_star_SJ_not_in_classinput_priority_Align.txt`: splice junctions in the STAR `SJ.out.tab` file that cannot be found in `class_input_priorityAlign.tsv`
* `in_star_SJ_not_in_classinput_priority_Chim.txt`: chimeric junctions in the STAR `SJ.out.tab` file that cannot be found in class_input_priorityChimeric.tsv

### Log files
There is a file called `wrapper.log` in the folder for every pipeline run, as well as for every sample. The goal of these files it to make it easier to look at the output from the jobs you submit with the pipeline by collecting it all in the same place. For example, the folder `output/GSE109774_colon_cSM_10_cJOM_10_aSJMN_0_cSRGM_0` will contain a `wrapper.log` file which has the `.out` and `.err` files concatenated for every job and every sample in that run; these outputs are sorted by job type (so the outputs for the mapping jobs for each sample will be next to each other, etc). There is also a `wrapper.log` file in each sample sub-folder; for example, `output/GSE109774_colon_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/SRR6546273` will contain this file. It contains the output for all `.out` and `.err` outputs from all the jobs run on that specific sample. The `wrapper.log` files are rewritten every time the pipeline is run on a sample.

## A note on annotation

## Decisions to make in the future

Should our pipeline be based on independent alignment of both reads, or alignment in the paired-end read mode? This will depend on which is better for detecting junctions.

## Processed Files annotation columns

For the processed files:

`splice_ann`: both sides of the junction are annotated as exon boundaries and the junction itself is found in the GTF; subset of `both_ann`

`both_ann`:  Both sides of the junction are annotated as exon boundaries; this is equivalent to `exon_annR1A` AND `exon_annR1B`

`just_both_ann`: Both sides of the junction are annotated as exon boundaries, but the splice is not annotated; `just_both_ann` UNION `splice_ann` = `both_ann`

`exon_annR1A`: The exon on the first half of the junction is at an annotated boundary

`exon_annR1B`: The exon on the second half of the junction is at an annotated boundary

`one_ann`: Exactly one of `exon_annR1A` and `exon_annR1B` is true

`none_ann`: Neither of the boundaries of the splice junction are annotated exon boundaries. 

`none_ann_known_gene`: `none_ann` is TRUE and `geneR1A_uniq` is NOT unknown or blank

`none_ann_unknown_gene`: `none_ann` is TRUE and `geneR1A_uniq` IS unknown or blank
