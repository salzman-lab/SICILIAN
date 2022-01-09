import pandas as pd

def get_gene_id(row):
#  return row["attribute"].split(";")[0].split()[1][1:-1]
  if "gene_name" in row["attribute"]:
    return row["attribute"].split("gene_name")[-1].split('"')[1]
  elif ";gene=" in row["attribute"]:
    return row["attribute"].split(";gene=")[-1].split(";")[0]
  if "gene_id" in row["attribute"]:
    return row["attribute"].split("gene_id")[-1].split('"')[1]

def round_down(num, divisor): 
    return num - (num%divisor)

# This is a class to create an annotator object that you can create based on a gtf file,
# that allows you to put in a chromosome and position and get all gene names in that area.

# Usage example:
# import sys 
# sys.path.insert(0, '/scratch/PI/horence/JuliaO/single_cell/scripts/') 
# import annotator
# 
# ann = annotator.Annotator(/scratch/PI/horence/JuliaO/single_cell/STAR_output/mm10_files/mm10.gtf) # this step can take a while - like 2 minutes
# ann.get_name_given_locus("chr1", 1003024) # returns all gene names separated by ","; if none, returns ""

class Annotator:
  def __init__(self, gtf_file, jump = 10000):
    self.jump = jump
    self.gtf_file = gtf_file
    self.unknown = "unknown"
    self.unknown_strand = "?"
    self.get_gtf_dict()

  def get_gtf_dict(self):
    # load in gtf
    gtf_df = pd.read_csv(self.gtf_file,sep="\t",names=["seqname","source","feature","start","end","score","strand","frame","attribute"],comment="#")
    # make gene id column
    gtf_df["gene_id"] = gtf_df.apply(get_gene_id, axis=1)
    
    # figure out how long to make each chromosome entry
    seqname_len_dict = {}
    for seqname in gtf_df["seqname"].unique():
        print(seqname)
        seqname_len_dict[seqname] = max(gtf_df[gtf_df["seqname"] == seqname]["end"])
        if seqname_len_dict[seqname] < max(gtf_df[gtf_df["seqname"] == seqname]["start"]):
            print("start more than end")
  
    # set up gtf dict to have a dictionary for each chromsome with entries for every "jump" in its length
    gtf_dict = {s : {r : {} for r in range(0, seqname_len_dict[s],self.jump)} for s in seqname_len_dict.keys()}
  
    # assign genes to their requisite ranges
    for seqname in seqname_len_dict:
        seqname_df = gtf_df[gtf_df["seqname"] == seqname]
        for gene_id in seqname_df["gene_id"].unique():
            if gene_id is not None:
              gene_df = seqname_df[seqname_df["gene_id"] == gene_id]
              if len(gene_df["strand"].unique()) == 1:
                strand = gene_df["strand"].unique()[0]
              else:
                strand = self.unknown_strand
    
              # assign gene to all ranges it falls within
              try:
                start = min(gene_df["start"])
              except:
                print("gene_id",gene_id)
                print("start failed", gene_df)
              try: 
                end = max(gene_df["end"])
              except:
                print("gene_id",gene_id)
                print("end failed",gene_df)
              for j in range(round_down(start,self.jump),round_down(end + self.jump, self.jump),self.jump):
                  gtf_dict[seqname][j][gene_id] = [start,end, strand] 
    self.gtf_dict = gtf_dict

  def get_name_given_locus(self, seqname, position, read_strand = "", stranded_library = False): 
   
      try: 
          poss_genes = self.gtf_dict[seqname][round_down(position,self.jump)] 
      except Exception as e: 

          if seqname not in self.gtf_dict.keys(): 
              if stranded_library:
                return self.unknown,  read_strand

              else:
                return self.unknown, self.unknown_strand 
          if position > max(self.gtf_dict[seqname].keys()): 
              if stranded_library:
                return self.unknown, read_strand
              else:
                return self.unknown, self.unknown_strand
          else: 
              raise e 
      if len(poss_genes) == 0: 
          if stranded_library:
              return self.unknown, read_strand 
          else:
              return self.unknown, self.unknown_strand 

      gene_names = [] 
      strands = []
      for gene, pos in poss_genes.items(): 
          if pos[0] <= position <= pos[1]:
              if stranded_library:
                  if pos[2] == read_strand:
                      gene_names.append(gene)
                      strands.append(pos[2])
              else:
                  gene_names.append(gene)
                  strands.append(pos[2])
      if len(gene_names) == 0:
          gene_names.append(self.unknown)

      if len(set(strands)) == 1:
          strand = strands[0]
      elif stranded_library:
          strand = read_strand
      else:
          strand = self.unknown_strand 
      return ",".join(gene_names), strand
