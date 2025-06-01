import os
import re
import sys
import subprocess
import pickle
import argparse
import itertools

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.dataset import Dataset
from utils.basicTypes import CompRes

def createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, outFilenamePattern, logFilename, 
						ucscName_other, alpha, overwriteFiles):

	python_module = "python/3.11.4"
	python_env    = "~/python-env-2024"

	jobName = f"alpha{alpha}.taus"
	
	# python3 5_sampleEvolTimes.py -sp_ucsc_name [UCSC name] -chr [chromosome name] -alpha [any number > 1] -cores [nb. of cores] [--overwrite, optional]
	# evolTimes-ests.mm39.chr16.alpha1.1.winSize3805.pickle

	with open(slurmFilename, "w") as slurm_file:
		codeLines  = f"module load {python_module}\nsource {python_env}/bin/activate\nstart_time=$(date +%s)\n\n"
		codeLines += f"python3 {dirCode}/5_sampleEvolTimes.py -sp_ucsc_name {ucscName_other} -alpha {alpha} -cores {cores} {'--overwrite' if (overwriteFiles) else ''}\n"
		codeLines += "\nend_time=$(date +%s)\ntotal_time=$(( end_time - start_time ))\n\ndeactivate\n"
		codeLines += 'echo "Job ID: $SLURM_JOB_ID"\n'
		codeLines += f"ls -l {outFilenamePattern}" 
		codeLines += " | awk '{sum += $5} END {printf \"Total disk usage: %.3f GB\\n\", sum / 1073741824}'\n"
		codeLines += 'echo "\nTotal time to run: $total_time seconds"\nsstat --allsteps -j $SLURM_JOB_ID --format=JobID,MaxRSS\n'
		slurm_file.write(f"#!/bin/bash\n#SBATCH -p compute\n#SBATCH -t {duration}\n#SBATCH --output={logFilename}\n#SBATCH --job-name={jobName}\n#SBATCH --mem={memory}G\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={cores}\n\n{codeLines}\n")
	return slurmFilename

####################################
# MAIN.
####################################
# Usage:   python3 5_sampleEvolTimes_runAll.py --cores [nb. of cores] [--overwrite, optional]
# Example: python3 ~/code/cluster/5_sampleEvolTimes_runAll.py --cores 80 --overwrite
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="[Cluster Version] Evolutionary Time Inference: Sampling.")
	parser.add_argument("--cores", help="Number of cores to be used during setup.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")
	
	args = parser.parse_args()
	overwriteFiles = args.overwrite
	nbcores = int(args.cores)

	# This data structure contains all needed information about the dataset.
	my_dataset  = Dataset()

	# Directory organization.
	dirCode  = "/home/p/priscila-do/code" 
	dirSlurm = "/home/p/priscila-do/code/slurm"

	dirWin   = my_dataset.dirWindows
	dirTemp  = my_dataset.dirTemp
	dirEstEvolTimes  = my_dataset.dirEstEvolTimes
	dirSampEvolTimes = my_dataset.dirSampEvolTimes

	print("****************************************************************")
	print("*            Sample evolutionary times - CLUSTER               *")
	print("* Process the posterior evolutionary times for the entire      *")
	print("* dataset. Sample evolutionary times for each window across    *") 
	print("* all chromosomes. The sampled evolutionary times are then     *")
	print("* used to derive distributions for evolutionary time and PCS   *")
	print("* size at both the whole-genome level and for each chromosome. *")
	print("****************************************************************")
	print(f"- Number of cores: {nbcores}")
	print(f"- Directory for Python code: {dirCode}")
	print(f"- Directory for Slurm scripts: {dirSlurm}")
	print(f"- Directory for windows (input): {dirWin}")
	print(f"- Directory for estimates (input): {dirEstEvolTimes}")
	print(f"- Directory for sampled info (output): {dirSampEvolTimes}")
	print(f"- Directory for temporary files: {dirTemp}")
	print(f"- Overwrite output files? {overwriteFiles}")

	# Check if all directories exist.
	for dirPath in [dirCode, dirSlurm, dirWin, dirEstEvolTimes, dirSampEvolTimes, dirTemp]:
		if (not os.path.exists(dirPath)):
			print(f"ERROR! Directory not found ({dirPath}).")
			sys.exit()		

	# Check if input files exist.
	winFilenames = [my_dataset.getOutFilename_computeWindows(qChrom) for qChrom in my_dataset.chromLst]
	abortRun     = (not all(os.path.isfile(filename) for filename in winFilenames))
	if (abortRun):
		print(f"ERROR! Check if window files exist for all chromosomes.")
		sys.exit()

	# Method information.
	alphas    = my_dataset.alphas
	alphas    = [1.1]
	UCSCnames = my_dataset.speciesUCSCnames

	# For each alpha and each species being compared to the reference genome.
	for (alpha, ucscName_other) in itertools.product(alphas,UCSCnames):

		# Slurm Parameters.
		comp_res = my_dataset.getCompRes_sampleEvolTimes(ucscName_other)
		duration = comp_res.time # hh:mm:ss
		memory	 = comp_res.mem  # GB
		
		logFilename = my_dataset.getLogFilename_sampleEvolTimes(ucscName_other,alpha)
		outFilenamePattern = my_dataset.getOutFilenamePattern_sampleEvolTimes(ucscName_other,alpha)

		# Create slurm script.
		print(f"[α={alpha}, {ucscName_other}] Create slurm script...")
		slurmFilename = os.path.join(dirSlurm,f"alpha{alpha}.{ucscName_other}.ests.slurm")
		createSlurmScript(	slurmFilename, memory, nbcores, duration, 
							dirCode, outFilenamePattern, logFilename, 
							ucscName_other, alpha, overwriteFiles)
		# Run slurm script.
		print(f"[α={alpha}] Running slurm script... {slurmFilename}")
		process = subprocess.Popen(['sbatch', slurmFilename]) 
