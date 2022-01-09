import annotator
import argparse
import pandas as pd
import pickle
import time
import os

def get_args():
  parser = argparse.ArgumentParser(description="create annotator for assembly")
  parser.add_argument("-g", "--gtf_path", help="the path to the gtf file to use for annotation", default=False)
  parser.add_argument("-a", "--assembly", help="The name of the assembly to pre-load annotation (so, mm10 for the 10th mouse assembly)")
  args = parser.parse_args()
  return args

def get_transcript_id(row):
    return row["attribute"].split("transcript_id")[-1].split('"')[1]

def get_exon_number(row):
  if "exon_number" in row["attribute"]:
    return int(row["attribute"].split("exon_number")[-1].split("\"")[1].split(";")[0])

#def get_exon_number2(row):
#  return int(row["attribute"].split("ID=exon")[1].split(";")[0].split("-")[-1])
#  if "gene_name" in row["attribute"]:
#    return int(row["attribute"].split("exon_number")[-1].split('"')[1])


def get_gtf(gtf_path):
  gtf_df = pd.read_csv(gtf_path,sep="\t",names=["seqname","source","feature","start","end","score","strand","frame","attribute"],comment="#")
  gtf_df = gtf_df[gtf_df["feature"] == "exon"]
  return gtf_df

def get_exon_bounds(gtf_df):
  exon_bounds = {}
  for name, group in gtf_df.groupby("seqname"):
    exon_bounds[name] = set()
    exon_bounds[name].update(set(group["start"]))
    exon_bounds[name].update(set(group["end"]))
  return exon_bounds

def get_splices(gtf_df):
  gtf_df["transcript_id"] = gtf_df.apply(get_transcript_id, axis=1)
#  try:
  gtf_df["exon_number"] = gtf_df.apply(get_exon_number, axis=1)

  splices = {}
  count = 0
  t0 = time.time()
  for name1, group1 in gtf_df[gtf_df["feature"] == "exon"].groupby("seqname"):
    splices[name1] = set()
    count += 1
  #  print(count, name1, time.time())
    for name2, group2 in group1.groupby("transcript_id"):
      for i in range(1,int(max(group2["exon_number"]))):
        splices[name1].add(tuple(sorted([group2[group2["exon_number"] == i].iloc[0]["end"],group2[group2["exon_number"] == i + 1].iloc[0]["start"]])))
  return splices

def main():

  args = get_args()
  gtf_df = get_gtf(args.gtf_path)

  exon_bounds = get_exon_bounds(gtf_df)
  pickle.dump(exon_bounds, open("{}/{}_exon_bounds.pkl".format(os.getcwd(), args.assembly), "wb"))
  print("dumped exon boundary annotator to {}/{}_exon_bounds.pkl".format(os.getcwd(), args.assembly))

  splices = get_splices(gtf_df)
  pickle.dump(splices, open("{}/{}_splices.pkl".format(os.getcwd(), args.assembly), "wb"))
  print("dumped splice junctions annotator to {}/{}_splices.pkl".format(os.getcwd(), args.assembly))
  
  ann = annotator.Annotator(args.gtf_path)
  pickle.dump(ann, open("{}/{}.pkl".format(os.getcwd(), args.assembly), "wb"))
  print("dumped annotator to {}/{}.pkl".format(os.getcwd(), args.assembly))
main()
