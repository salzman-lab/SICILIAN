# View Sequences
# Created: 11 September 2019
# Julia Olivieri
# Print sequences to the screen with colors as they are aligned

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import time

def print_seq_colored(string):
  bases = {"A" : 41, "T" : 43, "C" : 46, "G" : 44}
  for c in string:
    if c in bases:
      print("\033[1;30;{}m{}".format(bases[c],c), end="")
    else:
      print("\033[1;30;47m{}".format(c), end="")
  print()

def write_seq(juncName, out_file, seq_df, out = "True"):
  
  if out:
    out = open(out_file, "w")
    out.write(juncName + "\n")
  else:
    print(juncName)
  A_min = min(seq_df["posR1A"])
  B_min = min(seq_df["posR1B"])
  if len(set(seq_df["posR1A"])) > 1:
    for i,row in seq_df.iterrows():
      if out:
        out.write(" "*int(row["posR1A"] - A_min) + row["seqR1"][:1 + int(abs(row["posR1A"] - row["juncPosR1A"]))]+ "| |" + row["seqR1"][1 + int(abs(row["posR1A"] - row["juncPosR1A"])):] + "\n")
      else:
        string = " "*int(row["posR1A"] - A_min) + row["seqR1"][:1 + int(abs(row["posR1A"] - row["juncPosR1A"]))]+ "| |" + row["seqR1"][1 + int(abs(row["posR1A"] - row["juncPosR1A"])):] + " "
        print_seq_colored(string)
#         print(" "*int(row["posR1A"] - A_min) + row["seqR1"][:1 + int(abs(row["posR1A"] - row["juncPosR1A"]))]+ "| |" + row["seqR1"][1 + int(abs(row["posR1A"] - row["juncPosR1A"])):])
  elif len(set(seq_df["posR1B"])) > 1:
    for i,row in seq_df.iterrows():
      if out:
        out.write(" "*int(row["posR1B"] - B_min) + row["seqR1"][:1 + int(abs(row["posR1B"] - row["juncPosR1B"]))] + "| |" + row["seqR1"][1 + int(abs(row["posR1B"] - row["juncPosR1B"])):]+"\n")
      else:
        string = " "*int(row["posR1B"] - B_min) + row["seqR1"][:1 + int(abs(row["posR1B"] - row["juncPosR1B"]))] + "| |" + row["seqR1"][1 + int(abs(row["posR1B"] - row["juncPosR1B"])):] + " "
        print_seq_colored(string)
#         print(" "*int(row["posR1B"] - B_min) + row["seqR1"][:1 + int(abs(row["posR1B"] - row["juncPosR1B"]))] + "| |" + row["seqR1"][1 + int(abs(row["posR1B"] - row["juncPosR1B"])):])
  else:
    for i,row in seq_df.iterrows():
      if out:
        out.write(row["seqR1"] + "\n")
      else:
        string = row["seqR1"] + " " 
        print_seq_colored(string)
#         print(row["seqR1"])
  if out:
    out.write("\n" + "-"*100 + "\n\n")
    out.close()
  else:
#     print("\n" + "-"*100 + "\n")
    string = "\n" + "-"*100 + " \n"
    print_seq_colored(string)
  
def view_seq_CI(CI, juncName, num_seqs = 100):
  seq_df = pd.DataFrame(columns=["seqR1", "posR1A", "posR1B", "juncPosR1A", "juncPosR1B"])
  seq_df = seq_df.append(CI[CI["refName_newR1"] == juncName][seq_df.columns])
  try:
    seq_df = seq_df.sample(n=num_seqs)
  except:
    print("too few samples")
  seq_df = seq_df.sort_values(by=["posR1A","posR1B"])
  return seq_df

def view_seq(files, juncName):
  seq_df = pd.DataFrame(columns=["seqR1", "posR1A", "posR1B", "juncPosR1A", "juncPosR1B"])
  for i in range(len(files)):
    c = files[i] + "class_input_WithinBAM.tsv"
    CI = pd.read_csv(c, sep = "\t")
    seq_df = seq_df.append(CI[CI["refName_newR1"] == juncName][seq_df.columns])
  seq_df = seq_df.sort_values(by=["posR1A","posR1B"])
  return seq_df

def main():
  t0 = time.time()
  name = "TSP1_lung_1_S16_L003"
  usecols = ["seqR1", "posR1A", "posR1B", "juncPosR1A", "juncPosR1B"]
#  df = pd.read_csv("/scratch/PI/horence/Roozbeh/single_cell_project/output/TS_pilot_10X_withinbam_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/{}/class_input_WithinBAM.tsv".format(name), sep = "\t", 
#        usecols = ["refName_newR1", "seqR1", "posR1A", "juncPosR1A", "juncPosR1B", "posR1B"],
#        dtype = {"juncPosR1A" : "int32",
#                                "juncPosR1B" : "int32", "refName_newR1" : "category", "posR1A" : "int32"})
#  pickle.dump(df, open("/oak/stanford/groups/horence/JuliaO/pickled/tenX_CI_df.pkl", "wb"))
  df = pickle.load(open("/oak/stanford/groups/horence/JuliaO/pickled/tenX_CI_df.pkl", "rb"))

#  name = "B107921_M23_S267"
#  df = pd.read_csv("/scratch/PI/horence/Roozbeh/single_cell_project/output/TS_pilot_smartseq_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/{}/class_input_WithinBAM.tsv".format(name), sep = "\t")
  print("loaded df",name, time.time() - t0)
  juncName = "chr11:H19_3,Y_RNA,Metazoa_SRP,SNORA7,snoU13,SNORA1,SCGB1A1,CTD-2531D15.4:62422408:?|chr11:H19_3,Y_RNA,Metazoa_SRP,SNORA7,snoU13,SNORA1,SCGB1A1,CTD-2531D15.4:62423069:?|lin"
  while juncName != "done":
    juncName = input("juncName to view (write 'done' to stop):")
#    try:
    seq_df = view_seq_CI(df, juncName, num_seqs = 25)
    write_seq(juncName, "out.txt", seq_df, out = False)
#    except:
#      print("juncName not found")
  
main()
