import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import pandas as pd
import pyarrow
import time
import itertools
import math
from matplotlib import collections as matcoll
import numpy as np
import pickle
from statsmodels.stats.proportion import proportion_confint

def plot_conc(temp_df,outpath,dataset):
  print(dataset)
  temp_df["total_numReads"] = temp_df["refName_newR1"].map(temp_df.groupby("refName_newR1")["numReads"].sum())
  temp_df = temp_df.drop_duplicates(["refName_newR1","individual"])
  len_x = 100
  temp_df["bin"] = 0
  x_vals = np.geomspace(1, temp_df["total_numReads"].max() + 1, num=len_x)
  for i in range(len(x_vals) - 1):
#    print(i)
    temp_df.loc[(temp_df["total_numReads"] >= x_vals[i]) & (temp_df["total_numReads"] < x_vals[i + 1]),"bin"] = i + 1

  # neither, one, both
  neith_one_both = [[],[],[]]
  numReads = []
  for name, group in temp_df.groupby("bin"):
    numReads.append(x_vals[name - 1])
    for i in range(len(neith_one_both)):
      neith_one_both[i].append(0)
    for ref, group2 in group.groupby("refName_newR1"):
      neith_one_both[group2["inc_emp.p"].sum()][-1] += 1
#    for i in range(len(neith_one_both)):
#      neith_one_both[i][-1] = neith_one_both[i][-1]/group["refName_newR1"].nunique()

  label_dict = {0 : "failed in both", 1 : "passed in one, failed in one", 2 : "passed in both"}
  color_dict = {0 : "blue", 1 : "orange", 2 : "blue"}
  fill_dict = {0 : 'none',1 : "full", 2: "full"}
  for i in range(len(neith_one_both)):
    plt.plot(numReads,neith_one_both[i],marker="o",linestyle="",label = label_dict[i],color=color_dict[i],alpha=0.6,markersize=6,fillstyle=fill_dict[i])
    print(label_dict[i], sum(neith_one_both[i]))
    print(label_dict[i],"> 1", sum(neith_one_both[i][1:]))
  # plt.legend()
  plt.legend(bbox_to_anchor=(1, 1))
  
  plt.xscale("log")
  plt.yscale("log")

#  plt.ylim([-0.05,1.05])
  plt.xlabel("junction read count (binned)")
  plt.ylabel("fraction of junctions")
  savename = "{}conc_{}".format(outpath, dataset)
  plt.savefig(savename + ".png",bbox_inches = "tight")
  plt.savefig(savename + ".svg",format='svg',dpi=1200, bbox_inches = "tight")

  plt.close()

def get_dfs(individuals, dataset):
  ind_dfs = {}

  for ind in individuals:
    print(ind)
    df = pd.read_parquet("/scratch/groups/horence/JuliaO/single_cell/Process_Mouse_Lemur/scripts/output/Process_CI_10x/{}_{}_lung.pq".format(dataset, ind))
    df["individual"] = ind
    df["ind_numReads"] = df["refName_newR1"].map(df.groupby("refName_newR1")["numReads"].sum())
    df["10_cutoff"] = (df["ind_numReads"] >= 10)

    ind_dfs[ind] = df


  shared_juncs = set(ind_dfs[individuals[0]]["refName_newR1"]).intersection(set(ind_dfs[individuals[1]]["refName_newR1"]))
  df_10_1 = ind_dfs[individuals[0]]
  df_10_1 = df_10_1[df_10_1["ind_numReads"] > 1]
  df_10_2 = ind_dfs[individuals[1]]
  df_10_2 = df_10_2[df_10_2["ind_numReads"] > 1]
  shared_juncs_10 = set(df_10_1["refName_newR1"]).intersection(set(df_10_2["refName_newR1"]))
  both_passed = len(set(df_10_1[df_10_1["10_cutoff"]]["refName_newR1"]).intersection(set(df_10_2[df_10_2["10_cutoff"]]["refName_newR1"])))
  both_failed = len(set(df_10_1[~df_10_1["10_cutoff"]]["refName_newR1"]).intersection(set(df_10_2[~df_10_2["10_cutoff"]]["refName_newR1"])))

  print(dataset)
  print("number shared > 1", len(shared_juncs_10))
  print("number passed both 10 cutoff > 1",both_passed)
  print("number failed both 10 cutoff > 1",both_failed)
  print("frac consistent",(both_passed+both_failed)/len(shared_juncs_10))

  print(dataset)
  both_passed = len(set(ind_dfs[individuals[0]][ind_dfs[individuals[0]]["10_cutoff"]]["refName_newR1"]).intersection(set(ind_dfs[individuals[1]][ind_dfs[individuals[1]]["10_cutoff"]]["refName_newR1"])))
  both_failed = len(set(ind_dfs[individuals[0]][~ind_dfs[individuals[0]]["10_cutoff"]]["refName_newR1"]).intersection(set(ind_dfs[individuals[1]][~ind_dfs[individuals[1]]["10_cutoff"]]["refName_newR1"])))
  print("number shared:", len(shared_juncs))
  print("number passed both 10 cutoff:",both_passed)
  print("number failed both 10 cutoff:",both_failed)
  print("frac consistent",(both_passed+both_failed)/len(shared_juncs))

  df = ind_dfs[individuals[0]].append(ind_dfs[individuals[1]],ignore_index=True)
  df = df[df["refName_newR1"].isin(shared_juncs)]
  return df

def main():
  matplotlib.rcParams.update({'font.size': 16})
  outpath = "/scratch/PI/horence/JuliaO/single_cell/Methods_Paper/scripts/output/shared_concordance/"
  ind_dfs = {}
  individuals = ["Antoine","Stumpy"]
  dataset = "lemur4"
  df = get_dfs(individuals, dataset)
  plot_conc(df,outpath,dataset)

  individuals = ["P2","P3"]
  dataset = "HLCA4"
  df = get_dfs(individuals, dataset)
  plot_conc(df,outpath,dataset)


main()
