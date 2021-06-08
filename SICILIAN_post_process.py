# Wrapper script for STAR and statistical modeling
# Created by Julia Olivieri
# 17 June 2019

import glob
import os
import subprocess
import sys
import time
import argparse

def sbatch_file(file_name,out_path, name, job_name, time, mem, command, queue, dep="", dep_type = "afterok"):
  """Write sbatch script given parameters"""
  job_file = open(file_name, "w")
  job_file.write("#!/bin/bash\n#\n")
  job_file.write("#SBATCH --job-name=" + job_name + "\n")
  job_file.write("#SBATCH --output={}postprocess_log_files/{}.%j.out\n".format(out_path, job_name))
  job_file.write("#SBATCH --error={}postprocess_log_files/{}.%j.err\n".format(out_path, job_name))
  job_file.write("#SBATCH --time={}\n".format(time))
  #job_file.write("#SBATCH --qos=high_p\n")
  job_file.write("#SBATCH -p {}\n".format(queue))
#  job_file.write("#SBATCH --account=horence\n")
#  job_file.write("#SBATCH --partition=nih_s10\n")
  job_file.write("#SBATCH --nodes=1\n")
  job_file.write("#SBATCH --mem={}\n".format(mem)) 
  if dep != "":
    job_file.write("#SBATCH --dependency={}:{}\n".format(dep_type,dep))
    job_file.write("#SBATCH --kill-on-invalid-dep=yes\n")
  job_file.write("date\n")
  job_file.write(command + "\n")
  job_file.write("date\n")
  job_file.close()


def consolidate(out_path, run_name, data_format, exon_pickle_file, splice_pickle_file, queue, dep = ""):
  "Run the consolidate_GLM_output_files.R script to consolidated all GLM output files within an output directory into a single file"
  command = "Rscript scripts/consolidate_GLM_output_files.R {} {} {} {} ".format(out_path, run_name, exon_pickle_file, splice_pickle_file)
  if data_format == "10x":
    command += " 1 "
  else:
    command += " 0 "
  sbatch_file("run_consolidate.sh", out_path, run_name,"consolidate_{}".format(run_name), "48:00:00", "150Gb", command, queue, dep=dep)  # used 200Gb for CML 80Gb for others and 300 for 10x blood3
  return submit_job("run_consolidate.sh")

def process(out_path, run_name, gtf_file, exon_pickle_file, splice_pickle_file, queue, dep = ""):
  "Run the consolidate_GLM_output_files.R script to consolidated all GLM output files within an output directory into a single file"
  command = "python3 scripts/Process_CI_10x.py -d {} -o {} -g {} -e {} -s {}".format(out_path, run_name, gtf_file, exon_pickle_file, splice_pickle_file)
  sbatch_file("run_process.sh", out_path, run_name,"process_{}".format(run_name), "48:00:00", "350Gb", command, queue, dep=dep)  # used 200Gb for CML 80Gb for others and 300 for 10x blood3
  return submit_job("run_process.sh")

def postprocess(out_path, run_name, data_format, queue, dep = ""):
  "Run the consolidate_GLM_output_files.R script to consolidated all GLM output files within an output directory into a single file"
  command = "Rscript scripts/post_processing.R {} {} ".format(out_path, run_name)
  if data_format == "10x":
    command += " 1 "
  else:
    command += " 0 "

  sbatch_file("run_postprocess.sh", out_path, run_name,"postprocess_{}".format(run_name), "48:00:00", "200Gb", command, queue, dep=dep)  # used 200Gb for CML 80Gb for others and 300 for 10x blood3
  return submit_job("run_postprocess.sh")

def submit_job(file_name):
  """Submit sbatch job to cluster"""
  status, job_num = subprocess.getstatusoutput("sbatch {}".format(file_name))
  if status == 0:
    print("{} ({})".format(job_num, file_name))
    return job_num.split()[-1]
  else:
    print("Error submitting job {} {} {}".format(status, job_num, file_name))

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir', required=True, help='the output directory for the SICILIAN output files')
  parser.add_argument('-r', '--run', required=True, help='the name for the SICILIAN output files folder')
  parser.add_argument('-g', '--gtf', required=True, help='the path to the gtf file')
  parser.add_argument('-e', '--exon', required=True, help='the path to the exon pickle file')
  parser.add_argument('-s', '--splice', required=True, help='the path to the splice pickle file')
  parser.add_argument('-f', '--format', required=True, help='data format 10x or SS2', choices=["ss2","10x"])
  parser.add_argument('-q', '--queue', required=True, help='the queue for submitting jobs')
  parser.add_argument('-a', '--cons_step', required=True, help='the flag for running the consolidation step')
  parser.add_argument('-b', '--process_step', required=True, help='the flag for running the processing step')
  parser.add_argument('-c', '--postprocess_step', required=True, help='the flag for running the postprocessing model')
  args = parser.parse_args()


#################  Input arguments ###################
#  out_dir = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/output"
#  run_name = "HLCA_171205_tumor_10X"
#  gtf_file = "/oak/stanford/groups/horence/circularRNApipeline_Cluster/index/grch38_known_genes.gtf"
#  exon_pickle_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/scripts/STAR_wrapper/annotators/hg38_refseq_exon_bounds.pkl"
#  splice_pickle_file = "/oak/stanford/groups/horence/Roozbeh/single_cell_project/scripts/STAR_wrapper/annotators/hg38_refseq_splices.pkl"
#  data_format = "10x"
######################################################


#################  Input arguments ###################
  out_dir = args.dir
  run_name = args.run
  gtf_file = args.gtf
  exon_pickle_file = args.exon
  splice_pickle_file = args.splice
  data_format = args.format
  queue = args.queue
######################################################


######## flags for running postprocessing steps #########
  run_consolidate = args.cons_step
  run_process = args.process_step
  run_postprocess = args.postprocess_step
########################################################

  if data_format != "10x":
    run_process = False

  out_path = out_dir + "/{}/".format(run_name) 

  jobs = []
  job_nums = []

  if not os.path.exists("{}postprocess_log_files".format(out_path)):
    os.makedirs("{}postprocess_log_files".format(out_path))

  if run_consolidate:
    consolidate_jobid = consolidate(out_path, run_name, data_format, exon_pickle_file, splice_pickle_file, queue)
    jobs.append("consolidate_{}.{}".format(run_name, consolidate_jobid))
    job_nums.append(consolidate_jobid)
  else:
    consolidate_jobid = ""
  if run_process:
    process_jobid = process(out_path, run_name, gtf_file, exon_pickle_file, splice_pickle_file, queue, dep = ":".join(job_nums))
    jobs.append("process_{}.{}".format(run_name, process_jobid))
    job_nums.append(process_jobid)
  else:
    process_jobid = ""
  if run_postprocess:
    postprocess_jobid = postprocess(out_path, run_name, data_format, queue, dep = ":".join(job_nums))
    jobs.append("postprocess_{}.{}".format(run_name, postprocess_jobid))
    job_nums.append(postprocess_jobid)
  else:
    postprocess_jobid = ""
main()
