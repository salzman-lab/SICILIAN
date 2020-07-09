import numpy as np
import pandas as pd
import pyarrow

def get_tables(df, name, outpath):
  df["total_numReads"] = df["refName_newR1"].map(df.groupby("refName_newR1")["numReads"].sum())
  df["10_cutoff"] = (df["total_numReads"] >= 10)
  df["2_cutoff"] = (df["total_numReads"] >= 2)
  uniq_counts = df.groupby(["inc_emp.p","splice_ann","none_ann"])["refName_newR1"].nunique()
  uniq_10reads = df.groupby(["10_cutoff","splice_ann","none_ann"])["refName_newR1"].nunique()
  uniq_uniqreads = df.groupby(["uniq","splice_ann","none_ann"])["refName_newR1"].nunique()

  uniq_2reads = df.groupby(["2_cutoff","splice_ann","none_ann"])["refName_newR1"].nunique()
  uniq_table = {"Before filtering" : [], "After SICILIAN" : [], "After 10 read cutoff" : [], "After 2 read cutoff" : [], "After uniquely mapping cutoff" : []}
  uniq_table["Before filtering"].append(uniq_counts.loc[(True,True,False)] + uniq_counts.loc[(False,True,False)])
  uniq_table["Before filtering"].append(uniq_counts.loc[(True,False,True)] + uniq_counts.loc[(False,False,True)])
  uniq_table["After SICILIAN"].append(uniq_counts.loc[(True,True,False)])
  uniq_table["After SICILIAN"].append(uniq_counts.loc[(True,False,True)])
  uniq_table["After 10 read cutoff"].append(uniq_10reads.loc[(True,True,False)])
  uniq_table["After 10 read cutoff"].append(uniq_10reads.loc[(True,False,True)])
  uniq_table["After uniquely mapping cutoff"].append(uniq_uniqreads.loc[(True,True,False)])
  uniq_table["After uniquely mapping cutoff"].append(uniq_uniqreads.loc[(True,False,True)])

  uniq_table["After 2 read cutoff"].append(uniq_2reads.loc[(True,True,False)])
  uniq_table["After 2 read cutoff"].append(uniq_2reads.loc[(True,False,True)])
  uniq_df = pd.DataFrame.from_dict(uniq_table)
  uniq_df.set_index(np.array(["annotated","unannotated"]),inplace=True)
  for col in uniq_table.keys():
    if col not in ["Before filtering"]:
      uniq_df[col] = uniq_df[col].astype(str) + " (" + ((uniq_df[col]/uniq_df["Before filtering"]*100).round(decimals=2)).astype(str) + "%)"
  
  uniq_df.to_csv("{}{}_uniq.tsv".format(outpath, name),sep="\t")
  all_counts = df.groupby(["inc_emp.p","splice_ann","none_ann"])["numReads"].sum()
  all_10reads = df.groupby(["10_cutoff","splice_ann","none_ann"])["numReads"].sum()
  all_uniqreads = df.groupby(["uniq","splice_ann","none_ann"])["numReads"].sum()
#  all_uniq10reads = df.groupby(["uniq_10","splice_ann","none_ann"])["numReads"].sum()


  all_2reads = df.groupby(["2_cutoff","splice_ann","none_ann"])["numReads"].sum()
  all_table = {"Before filtering" : [], "After SICILIAN" : [], "After 10 read cutoff" : [], "After 2 read cutoff" : [], "After uniquely mapping cutoff" : []}
  all_table["Before filtering"].append(all_counts.loc[(True,True,False)] + all_counts.loc[(False,True,False)])
  all_table["Before filtering"].append(all_counts.loc[(True,False,True)] + all_counts.loc[(False,False,True)])
  all_table["After SICILIAN"].append(all_counts.loc[(True,True,False)])
  all_table["After SICILIAN"].append(all_counts.loc[(True,False,True)])
  all_table["After 10 read cutoff"].append(all_10reads.loc[(True,True,False)])
  all_table["After 10 read cutoff"].append(all_10reads.loc[(True,False,True)])
  all_table["After uniquely mapping cutoff"].append(all_uniqreads.loc[(True,True,False)])
  all_table["After uniquely mapping cutoff"].append(all_uniqreads.loc[(True,False,True)])

  all_table["After 2 read cutoff"].append(all_2reads.loc[(True,True,False)])
  all_table["After 2 read cutoff"].append(all_2reads.loc[(True,False,True)])
  all_df = pd.DataFrame.from_dict(all_table)
  all_df.set_index(np.array(["annotated","unannotated"]),inplace=True)
  for col in all_table.keys():
    if col not in ["Before filtering"]:
      all_df[col] = all_df[col].astype(str) + " (" + ((all_df[col]/all_df["Before filtering"]*100).round(decimals=2)).astype(str) + "%)"
  all_df.to_csv("{}{}_all.tsv".format(outpath, name),sep="\t")
  return uniq_df, all_df

def main():
  outpath = "/scratch/PI/horence/JuliaO/single_cell/Methods_Paper/scripts/output/count_tables/"

  datasets = {"HLCA4" : ["P2","P3"], "lemur4" : ["Stumpy","Antoine"]}
  ind_dict = {"P2" : "human 1","P3" : "human 2","Antoine" : "mouse lemur 1","Stumpy" : "mouse lemur 2"}

  # get HLCA tables
#  individuals = ["P1","P2","P3"]

  types = ["uniq","count"]
  out_dfs = {t :{"individual" : [], "# ann" : [], "# ann passing" : [], "frac ann passing" : [], "# unann" : [], "# unann passing" : [], "frac unann passing" : []} for t in types}

  fracs = {"uniq" : {"ann" : [], "unann" : []},"count" : {"ann" : [], "unann" : []}}
  for dataset, individuals in datasets.items():
    dfs = []
    for ind in individuals: 
      df = pd.read_parquet("/scratch/PI/horence/JuliaO/single_cell/Process_Mouse_Lemur/scripts/output/Process_CI_10x/{}_{}_lung.pq".format(dataset, ind))
#      dfs.append(df)
      uniq_df, all_df = get_tables(df, "{}_{}_lung".format(dataset, ind), outpath)
      
      type_dict = {"uniq" : uniq_df, "count" : all_df}
      for t in types:
        out_dfs[t]["individual"].append(ind_dict[ind])
        out_dfs[t]["# ann"].append(int(type_dict[t].loc["annotated","After 2 read cutoff"].split()[0]))
        out_dfs[t]["# ann passing"].append(int(type_dict[t].loc["annotated","After SICILIAN"].split()[0]))
        out_dfs[t]["frac ann passing"].append(int(type_dict[t].loc["annotated","After SICILIAN"].split()[0])/int(type_dict[t].loc["annotated","After 2 read cutoff"].split()[0]))
        out_dfs[t]["# unann"].append(int(type_dict[t].loc["unannotated","After 2 read cutoff"].split()[0]))
        out_dfs[t]["# unann passing"].append(int(type_dict[t].loc["unannotated","After SICILIAN"].split()[0]))
        out_dfs[t]["frac unann passing"].append(int(type_dict[t].loc["unannotated","After SICILIAN"].split()[0])/int(type_dict[t].loc["unannotated","After 2 read cutoff"].split()[0]))



      
      fracs["uniq"]["ann"].append(int(uniq_df.loc["annotated","After SICILIAN"].split()[0])/int(uniq_df.loc["annotated","After 2 read cutoff"].split()[0]))
      fracs["uniq"]["unann"].append(int(uniq_df.loc["unannotated","After SICILIAN"].split()[0])/int(uniq_df.loc["unannotated","After 2 read cutoff"].split()[0]))
      fracs["count"]["ann"].append(int(all_df.loc["annotated","After SICILIAN"].split()[0])/int(all_df.loc["annotated","After 2 read cutoff"].split()[0]))
      fracs["count"]["unann"].append(int(all_df.loc["unannotated","After SICILIAN"].split()[0])/int(all_df.loc["unannotated","After 2 read cutoff"].split()[0]))

  for v1 in ["uniq","count"]:
    for v2 in ["ann","unann"]:
      print(v1,v2,fracs[v1][v2],np.mean(fracs[v1][v2]))
  for t in types:
     pd.DataFrame(out_dfs[t]).to_csv("{}inc_{}.tsv".format(outpath,t),index=False,sep="\t")
     print("avg",t,"ann:",np.mean(out_dfs[t]["frac ann passing"]))
     print("avg",t,"unann:",np.mean(out_dfs[t]["frac unann passing"]))

   
#    all_df = pd.concat([x.reset_index(drop=True) for x in dfs], axis=1)
#    uniq_df, all_df = get_tables(all_df, dataset + "_full", outpath)

#  individuals = ["Bernard","Martine","Stumpy","Antoine"]
#  for ind in individuals:
#    df = pd.read_parquet("/scratch/PI/horence/JuliaO/single_cell/Process_Mouse_Lemur/scripts/output/Process_CI_10x/lemur_{}.pq".format(ind))
#    uniq_df, all_df = get_tables(df, "HLCA_" + ind, outpath)
#  df = pd.read_parquet("/scratch/PI/horence/JuliaO/single_cell/Process_Mouse_Lemur/scripts/output/Process_CI_10x/HLCA3.pq")
#  uniq_df, all_df = get_tables(df, "HLCA_full", outpath)

main()

