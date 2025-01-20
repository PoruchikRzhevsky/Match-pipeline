import os, sys, glob, re, yaml
import numpy as np
import pandas as pd

from process_data import *

sys.path.insert(1, ".")

procedures = ["coords_reproject", "matching", "filtering", "diagram",
			  "process"]

def run_all(sample, procedure, args, cluster=False):
	cwd = os.getcwd()
	if cluster:
		if procedure == "process":
			os.chdir(cwd)
			run(sample, "coords_reproject", args) ; os.chdir(cwd)
			run(sample, "matching", args) ; os.chdir(cwd)
			run(sample, "filtering", args) ; os.chdir(cwd)
			run(sample, "diagram", args) ; os.chdir(cwd)
		elif procedure not in procedures:
			print(f"{procedure.upper()} is not a procedure.")
			print("\nProcedures:\n" + " \n".join(procedures))
		else:
			run(sample, procedure, args)

def run(cluster, procedure, args):
	print(f"\nFor cluster {cluster} running procedure {procedure.upper()}.\n")
	
	# Set cluster directory
	obj_dir = f"{cluster}"
	
	# Load YAML file
	if not os.path.isdir(f'{obj_dir}'): os.system(f"mkdir -p {obj_dir}")
	if not os.path.isfile(f'{obj_dir}/{cluster}.yaml'): os.system(f"cp {cluster}.yaml {obj_dir}/")
	if not os.path.isdir(f'{obj_dir}/input'): os.system(f'mkdir -p {obj_dir}/input')
	if not os.path.isdir(f'{obj_dir}/plots'): os.system(f'mkdir -p {obj_dir}/plots')
	if not os.path.isdir(f'{obj_dir}/output'): os.system(f'mkdir -p {obj_dir}/output')
	if not os.path.isdir(f'{obj_dir}/output/coords'): os.system(f'mkdir -p {obj_dir}/output/coords')
	if not os.path.isdir(f'{obj_dir}/output/matched'): os.system(f'mkdir -p {obj_dir}/output/matched')
	os.chdir(obj_dir)
	y = yaml.load(open(f"{cluster}.yaml"), Loader=yaml.FullLoader)

	args = {"coords_reproject" : [y["cluster"], y["coords"], y["gaia_mag"], y["plots"], y["members"]],
			"matching" : [y["cluster"], y["gaia_mag"], y["matchrad"], y["trirad"], y["nobj"], y["plots"]],	
			"filtering" : [y["cluster"], y["plots"]],
			"diagram" : [y["cluster"], y["colour1"], y["colour2"], y["data_number"], y["adjust"], y["cmd"], y["test"], y["show_uncertainty"]]
			}

	eval(procedure)(*args[procedure])

if "__main__" == __name__:
	if len(sys.argv) < 3:
		print("Script for running procedures for cluster \n".upper())
		print("Usage: python3 run_all.py cluster/sample procedure")
		print("\nProcedures: " + " ".join(procedures))
		print("\nExample: python3 run_all.py NGC4649 spectra")
		print("Example: python3 run_all.py P21 reprocess")

	elif sys.argv[1]:
		cluster = sys.argv[1]
		procedure = sys.argv[2]
		args = sys.argv[3:]
		run_all(cluster, procedure, args, cluster=True)

	else:   
		print(f"Error: {sys.argv[1]} is invalid sample or cluster name")

