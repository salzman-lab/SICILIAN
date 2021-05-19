import argparse
from collections import defaultdict
import pandas as pd
import pyarrow
import pickle
import numpy as np
from glob import glob
import time
from tqdm import tqdm
from datetime import timedelta

def ensembl_name_map(gtf_file):
  gtf = pd.read_csv(gtf_file,sep="\t",header=None)
  gtf["ensembl_id"] = gtf[8].str.split("gene_id").str[1].str.split('"').str[1]
  gtf = gtf.drop_duplicates("ensembl_id")
  gtf["gene_name"] = gtf[8].str.split("gene_name").str[1].str.split('"').str[1]
  return gtf.set_index("ensembl_id")["gene_name"]

def add_genom_counts(df, lanes, data_path, ensembl_name):

  df["cell_gene_test"] = list(zip(df["geneR1A_uniq"], df["barcode"]))
  df["genom_gene_counts"] = np.nan

  gene_dfs = []
  for lane in lanes:
    try:
      print("{}{}/counts.tsv.gz".format(data_path, lane))
      gene_df = pd.read_csv("{}{}/counts.tsv.gz".format(data_path, lane),sep="\t")
      gene_df["gene_name"] = gene_df["gene"].map(ensembl_name)
      gene_dfs.append(gene_df)
    except Exception as e:
      print("raised exception",e)        
  gene_dfs = [d.groupby(["gene_name","cell"])["count"].sum() for d in gene_dfs]
  print("len gene dfs",len(gene_dfs))
  sum_df = sum(gene_dfs)
  for gd in gene_dfs:
    sum_df = sum_df.fillna(gd)
  print("sum_df",sum_df)
  
  try:
#    sum_df = sum_df.astype("Int32")
    idx = df[df["channel"] == lane[:-1]].index
    df.loc[idx,"genom_gene_counts"] = df.loc[idx]["cell_gene_test"].map(sum_df)
  except:
    print("no counts")
  return df


def get_args():
  parser = argparse.ArgumentParser(description="get method and number of files")
  parser.add_argument("-d","--data_paths",nargs="*",help="data paths for processing")
  parser.add_argument("-o","--outname",help="where to save output")
  parser.add_argument("-g","--gtf",help="gtf file for annotations")
  parser.add_argument("-e","--exon",help="pickle file for annotated exon boundaries")
  parser.add_argument("-s","--splice",help="pickle file for annotated splice junctions")
  parser.add_argument("--prefix",help="prefix for inc_emp.p",default="")
  parser.add_argument("--prefix2",help="prefix to help split many files",default="")
  parser.add_argument("--include_meta",action="store_true",help="include metadata")
  parser.add_argument("--suff",help="suffix for FDR",default="_fg_so_aa_ag_ae_il_0.15")


  args = parser.parse_args()
  return args


def get_names(data_paths, prefix):
  dp_dict = {}
  for data_path in data_paths:
    dp_dict[data_path] = {"err" : 0, "yes" : 0, "no" : 0, "names": []}
    for d in glob(data_path + "*/"):
      print(d)
      if prefix in d:
        try:
          with open(d + "class_input.tsv","r") as f:
            line = f.readline()[:-1].split("\t")
            # make sure file has the required column
            if "geneR1A_uniq" in line and "fileTypeR1" in line:
              dp_dict[data_path]["yes"] += 1
              dp_dict[data_path]["names"].append(d.split("/")[-2])
  
            else:
              dp_dict[data_path]["no"] += 1
        except Exception as e:
          print("exception", e)
          dp_dict[data_path]["err"] += 1
    print(len(dp_dict[data_path]["names"]),"names")
    print("dp_dict",dp_dict)
  return dp_dict

def process_lane(df,inc_refNames,meta_df,exon_bounds,splices,lane, include_meta):
#  if "fileTypeR1" in df.columns:
  df = df[df["fileTypeR1"] == "Aligned"]
  df = df.drop_duplicates(["barcode","UMI","refName_newR1"])

  df["barcode_refName"] = df["barcode"].astype(str) + df["refName_newR1"]
  barcode_name_vc = df["barcode_refName"].value_counts()
  name_vc = df["refName_newR1"].value_counts()
  df["numReads"] = df["barcode_refName"].map(barcode_name_vc)
  df = df.drop_duplicates(["refName_newR1","barcode"])
  df["gene_cell"] = df["geneR1A_uniq"].astype(str) + "_" + df["barcode"].astype(str)
#  df["gene_count_per_cell_no_filt"] = df["numReads"].groupby(df["gene_cell"]).transform("sum")
  df["gene_count_per_cell_no_filt"] = df["gene_cell"].map(df.groupby("gene_cell")["numReads"].sum())
  df["gene_frac_no_filt"] = df["numReads"]/df["gene_count_per_cell_no_filt"]
#  df["gene_frac_no_filt"] = df["numReads"]/df["gene_count_per_cell_no_filt"]
  df["gene_count_per_cell_filt"] = np.nan
  df["gene_frac_filt"] = np.nan
  for suffix in ["A","B"]:
      df["exon_annR1" + suffix] = False
      for name2, group in df.groupby("chrR1" + suffix):
        df.loc[group.index,"exon_annR1" + suffix] = group["juncPosR1" + suffix].isin(exon_bounds[name2])

  df["both_ann"] = (df["exon_annR1B"] & df["exon_annR1A"]).astype("bool")
  df["sort_junc"] = [tuple(sorted([x,y])) for x, y in zip(df.juncPosR1A,df.juncPosR1B)]
  df["splice_ann"] = False

  for name2, group in df.groupby("chrR1A"):
    sub_group = group[group["chrR1A"].astype(str) == group["chrR1B"].astype(str)]
    if name2 in splices:

      df.loc[sub_group.index,"splice_ann"] = sub_group["sort_junc"].isin(splices[name2])
  df["cell"] = lane + "_" + df["barcode"].astype(str)
  
  if include_meta:
    df = df.merge(meta_df, left_on="cell", right_on="cell", how = "left")
  df = df[[u for u in df.columns if u not in ["UMI","fileTypeR1","barcode_refName","gene_cell","sort_junc"]]]

  df["none_ann"] = ~df["splice_ann"] & ~df["both_ann"] & ~df["exon_annR1A"] & ~df["exon_annR1B"]
  df["one_ann"] = (df["exon_annR1A"] & ~df["exon_annR1B"]) | (~df["exon_annR1A"] & df["exon_annR1B"])
  df["just_both_ann"] = df["both_ann"] & ~df["splice_ann"]
  df["none_ann_known_gene"] = df["none_ann"] & (df["geneR1A_uniq"] != "unknown") & (~df["geneR1A_uniq"].isna())
  df["none_ann_unknown_gene"] = df["none_ann"] & ((df["geneR1A_uniq"] == "unknown") | (df["geneR1A_uniq"].isna()))

  print("end of processed lane")

  return df

def main():
  t0 = time.time()
  args = get_args()
  inc_refNames = []
  mult_lanes = False
  add_gen_count = False
  col = "emp.p_glmnet_constrained"
    
  gtf_file = args.gtf
  exon_bounds = pickle.load(open(args.exon,"rb"))
  splices = pickle.load(open(args.splice,"rb"))

  exon_bounds = defaultdict(set,exon_bounds)
  
  splices = defaultdict(set,splices)

  if add_gen_count:
    ensembl_name = ensembl_name_map(gtf_file)

  meta_df = [] 
  if args.include_meta:
    meta_df = pd.read_csv(args.data_paths[0] + "meta.tsv", sep="\t")

  use_cols = ["juncPosR1A","juncPosR1B","refName_newR1","geneR1A","UMI","barcode","fileTypeR1"]
  all_dfs = []
  dp_dict = get_names(args.data_paths, args.prefix2)
  print("dp_dict",dp_dict)

  for data_path in tqdm(args.data_paths):
    lane_dict = {}
    for name in tqdm(dp_dict[data_path]["names"]):
      lane = name[:-1]  #  name[:-1]
      print("name1",name)
      if lane not in lane_dict:
        lane_dict[lane] = []
      lane_dict[lane].append(name)
#  lane_dict = {'P1_3': ['P1_3_S1_L001', 'P1_3_S1_L002'], "P1_5" : ["P1_5_S3_L001"]}
    print("lane_dict",lane_dict)

    for lane in tqdm(lane_dict):

        dfs = []
        for name in lane_dict[lane]:
            print("name2",name, str(timedelta(seconds=time.time() - t0)))
            try:
              dfs.append(pd.read_csv("{}{}/class_input.tsv".format(data_path, name),
                                   dtype = {"p_predicted_glmnet" : "float16","p_predicted_glmnet_corrected" : "float16",
                                            "chrR1A" : "category", "chrR1B" : "category",
                                            "readClassR1" : "category", "geneR1A_uniq" : "category","geneR1A" : "category","geneR1B" : "category",
                                            "junc_cdf_glmnet" : "float16", "Organ" : "category",
                                           "splice_ann" : "bool", "both_ann" : "bool", "exon_annR1A" : "bool",
                                           "exon_annR1B" : bool},
                                   usecols = ["refName_newR1","UMI","barcode","geneR1A_uniq", 
                                              "juncPosR1A","juncPosR1B","chrR1A","chrR1B", "NHR1A","fileTypeR1"],
                                   sep = "\t"))
            except:
              print(name, "not processed")
    
        df = pd.concat(dfs)
        df.reset_index(drop=True,inplace=True)
        df = process_lane(df,inc_refNames,meta_df,exon_bounds,splices,lane, args.include_meta)
        print(df.columns)
        df["channel"] = lane
        if add_gen_count:
          df = add_genom_counts(df, lane_dict[lane], data_path, ensembl_name)
        print("processed lane")
        print("df.shape",df.shape)
        all_dfs.append(df)
 
  df = pd.concat(all_dfs)

  for c in df.columns:
    if str(df[c].dtype)[0] == "I":
      df[c] = df[c].astype("float32")
  df[[c for c in df.columns if c != "cell_gene_test"]].to_parquet(args.data_paths[0] + args.outname + args.prefix2+ ".pq")
  df[[c for c in df.columns if c != "cell_gene_test"]].to_csv(args.data_paths[0] + args.outname + args.prefix2+ ".tsv",sep="\t",index=False)
main()
