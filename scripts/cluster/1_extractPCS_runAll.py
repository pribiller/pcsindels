import os
import re
import sys
import subprocess
import pickle
import argparse

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.dataset import Dataset
from utils.basicTypes import CompRes

def createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, chainFilename, 
						outFilenames, logFilename, 
						dirTemp, ucscName_other):

	python_module = "python/3.11.4"
	python_env    = "~/python-env-2024"

	jobName = f"{ucscName_other}.pcs"

	chainFilenameTmp_compressed   = os.path.join(dirTemp, os.path.basename(chainFilename))
	chainFilenameTmp_decompressed = os.path.join(dirTemp, os.path.splitext(os.path.basename(chainFilename))[0])
	
	with open(slurmFilename, "w") as slurm_file:
		codeLines  = f"module load {python_module}\nsource {python_env}/bin/activate\nstart_time=$(date +%s)\n\n"
		codeLines += f"cp {chainFilename} {dirTemp}\n"
		codeLines += f"gzip -d {chainFilenameTmp_compressed}\n"
		codeLines += f"python3 {dirCode}/1_extractPCS.py -sp_ucsc_name {ucscName_other}\n"
		codeLines += f"rm {chainFilenameTmp_decompressed}\n"
		codeLines += "\nend_time=$(date +%s)\ntotal_time=$(( end_time - start_time ))\n\ndeactivate\n"
		codeLines += 'echo "Job ID: $SLURM_JOB_ID"\n'
		codeLines += f"ls -l {' '.join(outFilenames)}" 
		codeLines += " | awk '{sum += $5} END {printf \"Total disk usage: %.3f GB\", sum / 1073741824}'\n"
		codeLines += 'echo "\nTotal time to run: $total_time seconds"\nsstat --allsteps -j $SLURM_JOB_ID --format=JobID,MaxRSS\n'
		slurm_file.write(f"#!/bin/bash\n#SBATCH -p compute\n#SBATCH -t {duration}\n#SBATCH --output={logFilename}\n#SBATCH --job-name={jobName}\n#SBATCH --mem={memory}G\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={cores}\n\n{codeLines}\n")
	return slurmFilename

####################################
# MAIN.
####################################
# Usage:   python3 1_extractPCS_runAll.py
# Example: python3 ~/code/cluster/1_extractPCS_runAll.py --overwrite

if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="[Cluster Version] Extract PCSs (Perfectly Conserved Sequences) from Chains (ordered aligned blocks).")
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")
	
	args = parser.parse_args()
	overwriteFiles = args.overwrite

	# This data structure contains all needed information about the dataset.
	my_dataset  = Dataset()

	# Directory organization.
	dirCode     = "/home/p/priscila-do/code" 
	dirSlurm    = "/home/p/priscila-do/code/slurm"

	dirChains   = my_dataset.dirChains
	dirFasta    = my_dataset.dirFasta
	dirPCSs     = my_dataset.dirPCSs
	dirTemp     = my_dataset.dirTemp

	print("******************************************************")
	print("*   PCS (Perfectly conserved sequences) - CLUSTER    *")
	print("* Extract PCSs from chains (ordered aligned blocks). *")
	print("******************************************************")
	print(f"- Directory for Python code: {dirCode}")
	print(f"- Directory for Slurm scripts: {dirSlurm}")
	print(f"- Directory for chain files (input): {dirChains}") # e.g.: *.hg38.all.chain.gz
	print(f"- Directory for FASTA files (input): {dirFasta}")
	print(f"- Directory for PCSs (output): {dirPCSs}")
	print(f"- Directory for temporary files: {dirTemp}")
	print(f"- Overwrite output files? {overwriteFiles}")
	
	# Check if all directories exist.
	for dirPath in [dirCode, dirSlurm, dirChains, dirFasta, dirPCSs, dirTemp]:
		if (not os.path.exists(dirPath)):
			print(f"ERROR! Directory not found ({dirPath}).")
			sys.exit()		

	# Species information.
	ucscName_human   = my_dataset.refsp_ucscName
	speciesUCSCnames = my_dataset.speciesUCSCnames

	# For each species.
	for ucscName_other in speciesUCSCnames:

		# Check if input files exist.
		chainFilename = my_dataset.getChainFilename(ucscName_other,compressed=True)
		dirFasta_qry  = my_dataset.getDirFasta(ucscName_human)
		dirFasta_tar  = my_dataset.getDirFasta(ucscName_other)

		logFilename   = my_dataset.getLogFilename_extractPCS(ucscName_other)
		outFilenames  = [my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom) for qChrom in my_dataset.chromLst]

		skipRun = ((not overwriteFiles) and all(os.path.isfile(filename) for filename in outFilenames))
		if(skipRun): 
			print(f"WARNING! All files for species {ucscName_other} were already computed. Skipping run...")
			continue

		if ((chainFilename == None) or (dirFasta_qry == None) or (dirFasta_tar == None)):
			print(f"ERROR! {ucscName_other} : Check if chain file and FASTA files exist.")
			continue

		# Slurm Parameters.
		comp_res = my_dataset.getCompRes_extractPCS(ucscName_other)
		duration = comp_res.time # hh:mm:ss
		memory	 = comp_res.mem  # GB
		cores    = 1             # number of cores
		
		# Create slurm script.
		print(f"[{ucscName_other}] Create slurm script...")
		slurmFilename = os.path.join(dirSlurm,f"{ucscName_other}.pcs.slurm")
		createSlurmScript(	slurmFilename, memory, cores, duration, 
							dirCode, chainFilename, outFilenames, logFilename, dirTemp, ucscName_other)

		# Run slurm script.
		print(f"[{ucscName_other}] Running slurm script... {slurmFilename}")
		process = subprocess.Popen(['sbatch', slurmFilename])
		
