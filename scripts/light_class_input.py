import argparse
from collections import defaultdict
import datetime
import math
import numpy as np
import os
import pandas as pd
import pyarrow
import pickle
import pysam
import re
import time
import sys
#sys.path.insert(1, '/scratch/PI/horence/JuliaO/single_cell/STAR_wrapper/scripts/')
import annotator
from light_utils import *

def max_base(seq):
  base_counts = {"A" : [], "T" : [], "G" : [], "C" : []}
  for i in range(len(seq) - 9):
    for b in base_counts.keys():
      base_counts[b].append(seq[i:i + 10].count(b))
  return max(base_counts["A"]),max(base_counts["T"]),max(base_counts["G"]),max(base_counts["C"])

def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--bams', nargs="+",required=True, help='bams to parse (either one or two for paired end)')
  parser.add_argument("--outpath",help="folder to write output to")
  parser.add_argument("--UMI_bar", action="store_true", help="extract UMI and barcode")
  parser.add_argument("--paired", action="store_true", help="run once with each read primary and concatenate the files")
  parser.add_argument("--annotator", required=True, help="the path to the annotator pickle file")
  parser.add_argument("--gtf", required=True, help="the path to the gtf file")
  parser.add_argument("--stranded_library", action="store_true", help="Take gene strand information directly from the read strand")

  args = parser.parse_args()
  return args

def extract_info_align(CI_dict,bam_read,suffix,bam_file, ann, UMI_bar, stranded_library, fill_char = np.nan):
  sec_dict = {True: 0, False: 1}
  if UMI_bar:
    vals = bam_read.query_name.split("_")
    CI_dict["barcode"].append(vals[-2])
    CI_dict["UMI"].append(vals[-1])
  else:
    CI_dict["barcode"].append(fill_char)
    CI_dict["UMI"].append(fill_char)
  CI_dict["id"].append(bam_read.query_name)
  CI_dict["fileType" + suffix].append("Aligned")
  CI_dict["readLen" + suffix].append(bam_read.query_length)
  CI_dict["aScore{}A".format(suffix)].append(bam_read.get_tag("AS"))
  CI_dict["NH{}A".format(suffix)].append(bam_read.get_tag("NH"))
  CI_dict["HI{}A".format(suffix)].append(bam_read.get_tag("HI"))
  CI_dict["nmm{}A".format(suffix)].append(nmm(bam_read.get_tag("MD")))
  CI_dict["qual{}A".format(suffix)].append(bam_read.mapping_quality)
  cigar1, cigar2 = split_cigar_align(bam_read.cigarstring)
  M, S, I, D = get_SM(cigar1)
  CI_dict["M{}A".format(suffix)].append(M)
  CI_dict["S{}A".format(suffix)].append(S)
  M, S, I, D = get_SM(cigar2)
  CI_dict["M{}B".format(suffix)].append(M)
  CI_dict["S{}B".format(suffix)].append(S)
  CI_dict["cigar{}A".format(suffix)].append(cigar1)
  CI_dict["cigar{}B".format(suffix)].append(cigar2)
  seq = bam_read.query_sequence
  counts = count_stretch(seq)
  CI_dict["AT_run_" + suffix].append(max(counts["A"],counts["T"]))
  CI_dict["GC_run_" + suffix].append(max(counts["G"],counts["C"]))
  CI_dict["max_run_" + suffix].append(max(counts.values()))
  CI_dict["seq{}".format(suffix)].append(seq)
  CI_dict["entropy" + suffix].append(float('{:.3f}'.format(entropy(seq))))
  maxA, maxT, maxG, maxC = max_base(seq)
  CI_dict["maxA_10mer" + suffix].append(maxA)
  CI_dict["maxT_10mer" + suffix].append(maxT)
  CI_dict["maxG_10mer" + suffix].append(maxG)
  CI_dict["maxC_10mer" + suffix].append(maxC)

  refName, chrA, geneA, posA, chrB, geneB, posB = readObj_refname(bam_read.flag, bam_read.cigarstring, bam_file.get_reference_name(bam_read.tid), bam_read.reference_start + 1, ann, fill_char, stranded_library)
#   print("refName",refName)
  CI_dict["refName_AB" + suffix].append(refName)
  CI_dict["chr{}A".format(suffix)].append(chrA)
  CI_dict["chr{}B".format(suffix)].append(chrB)
  CI_dict["gene{}A".format(suffix)].append(geneA)
  CI_dict["gene{}B".format(suffix)].append(geneB)
  CI_dict["juncPos{}A".format(suffix)].append(int(posA))
  if np.isnan(posB):
    CI_dict["juncPos{}B".format(suffix)].append(posB)
  else:
    CI_dict["juncPos{}B".format(suffix)].append(int(posB))
  CI_dict["read_strand{}A".format(suffix)].append(read_strand(bam_read.flag))
  CI_dict["read_strand{}B".format(suffix)].append(fill_char)
  CI_dict["flag{}A".format(suffix)].append(bam_read.flag)
  CI_dict["flag{}B".format(suffix)].append(fill_char)

  CI_dict["primary{}A".format(suffix)].append(sec_dict[bam_read.is_secondary])
  CI_dict["primary{}B".format(suffix)].append(fill_char)
  empty_cols = ["aScore{}B".format(suffix),"qual{}B".format(suffix),"NH{}B".format(suffix),"nmm{}B".format(suffix), "HI{}B".format(suffix)]
  for c in empty_cols:
    CI_dict[c].append(fill_char)
  return CI_dict

def extract_info_chim(CI_dict,bam_read1,bam_read2,suffix, bam_file, ann, UMI_bar, stranded_library, fill_char = np.nan):
  assert bam_read1.query_name == bam_read2.query_name
  sec_dict = {True: 0, False: 1}
  if UMI_bar:
    vals = bam_read1.query_name.split("_")
    CI_dict["barcode"].append(vals[-2])
    CI_dict["UMI"].append(vals[-1])
  else:
    CI_dict["barcode"].append(fill_char)
    CI_dict["UMI"].append(fill_char)
  reads = [bam_read1,bam_read2]
  CI_dict["fileType" + suffix].append("Chimeric")
  halves = ["A","B"]
  CI_dict["id"].append(bam_read1.query_name)
  
  CI_dict["readLen" + suffix].append(bam_read1.query_length)
  
  seq = bam_read1.query_sequence
  counts = count_stretch(seq)
  CI_dict["AT_run_" + suffix].append(max(counts["A"],counts["T"]))
  CI_dict["GC_run_" + suffix].append(max(counts["G"],counts["C"]))
  CI_dict["max_run_" + suffix].append(max(counts.values()))
  CI_dict["seq{}".format(suffix)].append(seq)
  CI_dict["entropy" + suffix].append(float('%3.f'%(entropy(seq))))
  maxA, maxT, maxG, maxC = max_base(seq)
  CI_dict["maxA_10mer" + suffix].append(maxA)
  CI_dict["maxT_10mer" + suffix].append(maxT)
  CI_dict["maxG_10mer" + suffix].append(maxG)
  CI_dict["maxC_10mer" + suffix].append(maxC)

  refName, chrA, geneA, posA, chrB, geneB, posB  = chim_refName([x.flag for x in reads], [x.cigarstring for x in reads], [x.reference_start + 1 for x in reads], [bam_file.get_reference_name(x.tid) for x in reads], ann, stranded_library)
  CI_dict["refName_AB" + suffix].append(refName)
#   split_ref = refName.split("|")
  CI_dict["chr{}A".format(suffix)].append(chrA)
  CI_dict["chr{}B".format(suffix)].append(chrB)
  CI_dict["gene{}A".format(suffix)].append(geneA)
  CI_dict["gene{}B".format(suffix)].append(geneB)
  CI_dict["juncPos{}A".format(suffix)].append(int(posA))
  CI_dict["juncPos{}B".format(suffix)].append(int(posB))
  for i in range(2):
    CI_dict["aScore{}{}".format(suffix,halves[i])].append(reads[i].get_tag("AS"))
    CI_dict["qual{}{}".format(suffix,halves[i])].append(reads[i].mapping_quality)
    CI_dict["NH{}{}".format(suffix,halves[i])].append(reads[i].get_tag("NH"))
    CI_dict["HI{}{}".format(suffix,halves[i])].append(reads[i].get_tag("HI"))
    CI_dict["nmm{}{}".format(suffix,halves[i])].append(nmm(reads[i].get_tag("MD")))
    cigar = split_cigar_chim(reads[i].cigarstring)
    M, S, I, D = get_SM(cigar)
    CI_dict["M{}{}".format(suffix,halves[i])].append(M)
    CI_dict["S{}{}".format(suffix,halves[i])].append(S)
    CI_dict["cigar{}{}".format(suffix,halves[i])].append(cigar)
    CI_dict["read_strand{}{}".format(suffix,halves[i])].append(read_strand(reads[i].flag))
    CI_dict["flag{}{}".format(suffix,halves[i])].append(reads[i].flag)

    CI_dict["primary{}{}".format(suffix,halves[i])].append(sec_dict[reads[i].is_secondary])
  return CI_dict

def get_final_df(bam_files,j,suffixes,ann,UMI_bar,t0,gtf, stranded_library):
  CI_dfs = []
  for i in range(len(bam_files)):
    if i == 1:
      read_ids = set(CI_dfs[0]["id"])
    else:
      read_ids = set()
    suffix = suffixes[i]
    col_bases = ["aScore","M","S","nmm","qual","NH","HI","cigar", "juncPos", "gene", "chr", "read_strand","primary","flag"]
    columns = ["id","readLen" + suffix, "fileType" + suffix,"seq" + suffix,"AT_run_" + suffix,"GC_run_" + suffix,
               "max_run_" + suffix,"entropy" + suffix,"refName_AB" + suffix, "UMI","barcode","maxA_10mer" + suffix,"maxT_10mer" + suffix,"maxG_10mer" + suffix,"maxC_10mer" + suffix]
    for c in col_bases:
  #     for r in ["R1"]:
      for l in ["A","B"]:
        columns.append("{}{}{}".format(c,suffix,l))
    CI_dict = {c : [] for c in columns}
    count = 0
    first = False
    if i == 0:
      genomic_alignments = {}
    alignFile = pysam.AlignmentFile(bam_files[i])
    # columns
    for bam_read in alignFile.fetch(until_eof=True):
  #     suffix = "R1"

      # make sure read is mapped
      if not bam_read.is_unmapped:
        if (i == 0) or (not bam_read.is_secondary and bam_read.query_name in read_ids):
          # it's a chimeric alignment and we need another line from it
          if bam_read.has_tag("ch") and not first:
                prev_read = bam_read
                first = True
          else:

            # add info from chimeric read
            if bam_read.has_tag("ch"):
              count += 1

              # note: removing chim for this test ONLY; uncomment after
              CI_dict = extract_info_chim(CI_dict,prev_read,bam_read,suffix, alignFile, ann, UMI_bar, stranded_library)
              first = False

            # add info from align read
            elif "N" in bam_read.cigarstring:
              count += 1
              CI_dict = extract_info_align(CI_dict,bam_read,suffix,alignFile, ann, UMI_bar, stranded_library)

            # save genomic alignment information
            else:
              if i == 0:
                if bam_read.query_name not in genomic_alignments:
                  genomic_alignments[bam_read.query_name] = bam_read.get_tag("AS")
                else:
                  genomic_alignments[bam_read.query_name] = max(bam_read.get_tag("AS"), genomic_alignments[bam_read.query_name])
              else:
                CI_dict = extract_info_align(CI_dict,bam_read,suffix,alignFile, ann, UMI_bar, stranded_library)
  #     if count == 10000:
  #       continue
    CI_df = pd.DataFrame.from_dict(CI_dict)
    if i == 0:
      genomic_alignments = defaultdict(lambda: np.nan,genomic_alignments)

      CI_df["genomicAlignmentR1"] = 0
      # print(CI_df.loc[CI_df["id"].isin(genomic_alignments)].shape)
      CI_df.loc[CI_df["id"].isin(genomic_alignments),"genomicAlignmentR1"] = 1
      CI_df["genomic_aScoreR1"] = CI_df["id"].map(genomic_alignments)
      CI_df["spliceDist"] = abs(CI_df["juncPosR1A"] - CI_df["juncPosR1B"])
    CI_dfs.append(CI_df)
  if len(bam_files) == 2:
    final_df = pd.merge(left=CI_dfs[0],right=CI_dfs[1][[c for c in CI_dfs[1].columns if c not in ["UMI","barcode","seqR2"]]],how="left",left_on="id",right_on="id")
    final_df["read_strand_compatible"] = 1
    final_df.loc[final_df["read_strandR1A"] == final_df["read_strandR2A"],"read_strand_compatible"] = 0
    final_df["location_compatible"] = final_df.apply(get_loc_flag,axis=1)
  else:
    final_df = CI_dfs[0]
  #  final_df.fillna(np.nan,inplace=True)
  float_cols = ["aScoreR1A","nmmR1A","qualR1A","NHR1A","primaryR1A","genomic_aScoreR1","HIR1A"]
  if len(bam_files) == 2:
    float_cols += ["readLenR2","AT_run_R2",
                "GC_run_R2","max_run_R2","aScoreR2A","aScoreR2B","MR2A","MR2B","SR2A","SR2B","nmmR2A","nmmR2B",
                "qualR2A","qualR2B","NHR2A","NHR2B","juncPosR2A","juncPosR2B","primaryR2A","primaryR2B",  "HIR2A", "HIR2B"]

  #  for c in float_cols:
  #    final_df[c] = final_df[c].astype("Int32")
  print("started modify", time.time() - t0)
  #  for c in final_df.columns:
  #    str_dtype = final_df[c].dtype 
  #    if str(str_dtype)[0] == "i":
  #      final_df[c] = final_df[c].astype("I" + str_dtype[1:])

  final_df = modify_refnames(final_df, gtf, stranded_library) 

  print("ended modify", time.time() - t0)
  final_df["max_id_priority"] = final_df["id"].map(final_df.groupby("id")["HIR1A"].min())

  #  for c in final_df.columns:
  #    if str(final_df[c].dtype)[0] == "I":
  #      final_df[c] = final_df[c].astype("float32")

  #  final_df.to_parquet(args.outpath + "class_input_final.pq")
  final_df["primary_bam"] = j
  final_df = final_df[final_df["HIR1A"] == final_df["max_id_priority"]]
#  secondary = final_df[final_df["HIR1A"] != final_df["max_id_priority"]]

  return final_df

def main():
  t0 = time.time()
  args = get_args()

  bam_files = args.bams
  gtf = args.gtf
  print("bam_files",bam_files)

  annotator_path = args.annotator
  ann = pickle.load(open(annotator_path, "rb"))

  suffixes = ["R1","R2"]

  final_dfs = []
#  final_dfs_secondary = []
  
  if args.paired:
    n_rounds = 2
  else:
    n_rounds = 1

  for j in range(n_rounds):
    if j == 1:
      bam_files.reverse()
    print("bam_files",bam_files)
    primary = get_final_df(bam_files,j,suffixes,ann,args.UMI_bar,t0,gtf,args.stranded_library)
    final_dfs.append(primary)
#    final_dfs_secondary.append(secondary)
#  final_df = pd.concat(final_dfs,axis=0).reset_index(drop=True)


#  for c in final_df.columns:
#    if str(final_df[c].dtype)[0] == "I":
#      final_df[c] = final_df[c].astype("float32")

#  final_df.to_parquet(args.outpath + "class_input_final.pq")

  pd.concat(final_dfs,axis=0).reset_index(drop=True).to_parquet(args.outpath + "class_input.pq")
#  pd.concat(final_dfs_secondary,axis=0).reset_index(drop=True).to_parquet(args.outpath + "class_input_secondary.pq")

#  final_df[final_df["HIR1A"] == final_df["max_id_priority"]].to_hdf(args.outpath + "class_input.h5", key="class_input")
#  final_df[final_df["HIR1A"] != final_df["max_id_priority"]].to_hdf(args.outpath + "class_input_secondary.h5", key="class_input_secondary")

  pd.concat(final_dfs,axis=0).reset_index(drop=True).to_csv(args.outpath + "class_input.tsv",sep="\t",index=False)
#  pd.concat(final_dfs_secondary,axis=0).reset_index(drop=True).to_csv(args.outpath + "class_input_secondary.tsv",sep="\t",index=False)
#  final_df["juncPosR1A"] = final_df["juncPosR1A"].astype("Int32")
#  final_df["juncPosR1B"] = final_df["juncPosR1B"].astype("Int32") 

  print("total time",datetime.timedelta(seconds=time.time() - t0))
main()
