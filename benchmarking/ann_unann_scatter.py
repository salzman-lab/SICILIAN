import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import pandas as pd
import pyarrow
import time

def get_fracs(df, col, categories):
#   categories = ["splice_ann","both_ann","one_ann","none_ann"]
  count_out = {}
  uniq_out = {}
  ind_out = {"count" : {}, "uniq" : {}}
  
  for cat in categories:
    ind_out["uniq"][cat] = df[df[cat] & df[col]]["refName_newR1"].nunique()/df[df[cat]]["refName_newR1"].nunique()
    ind_out["count"][cat] = df[df[cat] & df[col]]["numReads"].sum()/df[df[cat]]["numReads"].sum()
  return ind_out

def get_totals(df, total_dict, categories):
  for cat in categories:
    total_dict["uniq"][cat].append(df[df[cat]]["refName_newR1"].nunique())
    total_dict["count"][cat].append(df[df[cat]]["numReads"].sum())
  return total_dict

def plot_totals(total_dict,categories,individuals,outpath):
  name_dict = {"Antoine" : "lemur 1", "Stumpy" : "lemur 2", "P2" : "human 1", "P3" : "human 2"}
  idx = range(len(name_dict))
  
  space = 0.8/len(categories)
  fill_dict = {"splice_ann" : True, "just_both_ann" : False, "one_ann" : False, "none_ann" : True}
  edge_color_dict = {"splice_ann" : 'none',"just_both_ann" : u'#ff7f0e',"one_ann" : u'#2ca02c',"none_ann" : 'none'}
  alpha_dict = {"splice_ann" : 0.6, "just_both_ann" : 1, "one_ann" : 1, "none_ann" : 0.6}
  yax_dict = {"uniq" : "number of junctions\nper category before filtering", "count" : "number of junctional reads\nper category before filtering"}
  cat_dict = {"splice_ann" : "junction annotated", "just_both_ann" : "both exons annotated (not splice)", "one_ann" : "one exon annotated", "none_ann" : "neither exon annotated", "both_ann" : "both exons annotated"}

  for y in ["uniq","count"]:
    for j in range(len(categories)):
      cat = categories[j]
      plt.bar([j*space + i - 0.3 for i in idx], total_dict[y][cat],width=space - 0.03, fill=fill_dict[cat],edgecolor=edge_color_dict[cat],linewidth=2,label=cat_dict[cat],alpha=0.6)
    plt.xticks(idx,[name_dict[ind] for ind in individuals])

    plt.ylabel(yax_dict[y])
    savename = "{}{}_total".format(outpath,y)
    plt.savefig(savename + ".png",bbox_inches = "tight")
    plt.savefig(savename + ".svg",format='svg',dpi=1200, bbox_inches = "tight")
    plt.legend(bbox_to_anchor=(2, 2))
    plt.savefig(savename + "_legend.svg",format='svg',dpi=1200, bbox_inches = "tight")

    plt.close() 


def plot_fracs(categories, types, cols, out_lists, individuals, outpath):
  color_dict = {"splice_ann" : u'#1f77b4',"just_both_ann" : u'#1f77b4',"one_ann" : u'#ff7f0e',"none_ann" : u'#ff7f0e'}
  edge_color_dict = {"splice_ann" : 'none',"just_both_ann" : u'#1f77b4',"one_ann" : u'#ff7f0e',"none_ann" : 'none'}

  fill_dict = {"splice_ann" : "full","just_both_ann" : "none","one_ann" : "none","none_ann" : "full"}
  name_dict = {"Antoine" : "lemur 1", "Stumpy" : "lemur 2", "P2" : "human 1", "P3" : "human 2"}
  y_dict = {"uniq" : "fraction of junctions passing\n", "count" : "fraction of junctional reads passing\n"}
  cat_dict = {"splice_ann" : "junction annotated", "just_both_ann" : "both exons annotated (not splice)", "one_ann" : "one exon annotated", "none_ann" : "neither exon annotated", "both_ann" : "both exons annotated"}
  col_dict = {"inc_emp.p" : "SICILIAN", "hard" : "hard filters", "uniq" : "uniquely mapping read filter", "10_cutoff" : "10 read filter", "uniq_10" : "uniquely mapping and 10 read filters"}
  alpha = 0.6
  xvals = range(len(categories))
  for col in cols:
    for t in types:
      for cat in categories:
#        plt.plot(xvals, out_lists[col][t][cat],alpha = alpha, label = cat_dict[cat],linestyle="",marker="o",markersize=8, color = color_dict[cat], fillstyle=fill_dict[cat],markeredgecolor=edge_color_dict[cat])
        plt.plot(xvals, out_lists[col][t][cat],alpha = alpha, label = cat_dict[cat],linestyle="",marker="o",markersize=8, fillstyle=fill_dict[cat])

      plt.xticks(xvals, [name_dict[i] for i in individuals],rotation="vertical")
      plt.ylabel(y_dict[t] + " " + col_dict[col] + " per category")
      fig = plt.gcf()
      fig.set_size_inches(2,4)
      plt.xlim([-0.5,len(categories) - 0.5])
      plt.ylim([-0.05,1.05])
      savename = "{}{}_{}".format(outpath, col, t)

      plt.savefig(savename + ".png",bbox_inches = "tight")
      plt.savefig(savename + ".svg",format='svg',dpi=1200, bbox_inches = "tight")
      plt.legend(bbox_to_anchor=(2, 2))
      plt.savefig(savename + "_legend.svg",format='svg',dpi=1200, bbox_inches = "tight")

      plt.close() 
  

def main():
  matplotlib.rcParams.update({'font.size': 16})

  outpath = "/scratch/PI/horence/JuliaO/single_cell/Methods_Paper/scripts/output/ann_unann_scatter/"
  t0 = time.time()
  types = ["count","uniq"]
  ind_prefix = {"Antoine" : "MLCA_ANTOINE_LUNG", "Stumpy" : "10X_P8_0", "P2" : "P2", "P3" : "P3"}
  data_loc = {"Antoine" : "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Lemur_Antoine_10X_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/", "Stumpy" : "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/Lemur_Stumpy_10X_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/", "P2" : "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/HLCA_180607_10X_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/", "P3" : "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output/HLCA_180607_10X_cSM_10_cJOM_10_aSJMN_0_cSRGM_0/"}

  individuals = ["P2","P3","Antoine","Stumpy"]
  categories = ["splice_ann","just_both_ann","one_ann","none_ann"]
  datasets = {"P2" : "HLCA4", "P3" : "HLCA4", "Antoine" : "lemur4", "Stumpy" : "lemur4"}
  cols = ["inc_emp.p","10_cutoff","uniq","uniq_10","hard"]
  total_dict = {k : {c : [] for c in categories} for k in ["uniq","count"]}
  out_lists = {x : {t : {c : [] for c in categories} for t in types} for x in cols}
  for ind in individuals:
    print(ind, time.time() - t0)
    df = pd.read_parquet("/scratch/groups/horence/JuliaO/single_cell/Process_Mouse_Lemur/scripts/output/Process_CI_10x/{}_{}_lung.pq".format(datasets[ind], ind))
    uniq_juncs = list(pd.read_csv("/scratch/PI/horence/JuliaO/single_cell/Methods_Paper/scripts/output/parse_uniquely_mapping/{}_all_counts.tsv".format(ind_prefix[ind]),names=["junctions","counts"],sep="\t")["junctions"])
    df["uniq"] = False
    df.loc[df["refName_newR1"].isin(uniq_juncs),"uniq"] = True
    inc_file = "{}{}_hard_inc_passed_fg_so_aa_ag_ae_il_1.txt".format(data_loc[ind], ind_prefix[ind])
    inc_refNames = set(list(pd.read_csv(inc_file,sep="\t",header=None,names=["refName_newR1"])["refName_newR1"]))
    df["hard"] = False
    df.loc[df["refName_newR1"].isin(inc_refNames),"hard"] = True

    df["total_numReads"] = df["refName_newR1"].map(df.groupby("refName_newR1")["numReads"].sum())
    print("total numreads", time.time() - t0)

    df = df[df["total_numReads"] >= 2]
    print("2 cutoff", time.time() - t0)

    df["10_cutoff"] = (df["total_numReads"] >= 10)
    print("10 cutoff", time.time() - t0)
    df["uniq_10"] = (df["10_cutoff"] & df["uniq"])
    total_dict = get_totals(df, total_dict, categories)
    print("\tloaded")
    for col in cols:
      ind_out = get_fracs(df, col, categories)
      for cat in categories:
        for t in types:
          out_lists[col][t][cat].append(ind_out[t][cat])

  plot_fracs(categories, types, cols, out_lists, individuals, outpath)

  plot_totals(total_dict,categories,individuals,outpath)
  print("total_dict",total_dict)
main()
