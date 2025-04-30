import os
# import argparse
import re
import sys
from datetime import datetime
import subprocess
import pickle

def createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, dirFasta_qry, dirFasta_tar, chainFilename, dirPCSs,
						ucscName_human, ucscName_other, minPCSsize):

	python_module = "python/3.11.4"
	python_env    = "~/python-env-2024"

	jobName = f"{ucscName_other}.pcs"

	chainFilenameTmp_compressed   = os.path.join(dirPCSs, os.path.basename(chainFilename))
	chainFilenameTmp_decompressed = os.path.join(dirPCSs, os.path.splitext(os.path.basename(chainFilename))[0])
	
	with open(slurmFilename, "w") as slurm_file:
		codeLines  = f"module load {python_module}\nsource {python_env}/bin/activate\nstart_time=$(date +%s)\n\n"
		codeLines += f"cp {chainFilename} {dirPCSs}\n"
		codeLines += f"gzip -d {chainFilenameTmp_compressed}\n"
		codeLines += f"python3 {dirCode}/1-extractPCS.py -refsp_ucscname {ucscName_human} -othsp_ucscname {ucscName_other} -chain_file {chainFilenameTmp_decompressed} -refsp_fastadir {dirFasta_qry} -othsp_fastadir {dirFasta_tar} -pcs_dir {dirPCSs} -pcs_minsize {minPCSsize}\n"
		codeLines += f"rm {chainFilenameTmp_decompressed}\n"
		codeLines += "\nend_time=$(date +%s)\ntotal_time=$(( end_time - start_time ))\n\ndeactivate\n"
		codeLines += 'echo "Job ID: $SLURM_JOB_ID"\necho "Total time to run: $total_time seconds"\n'
		slurm_file.write(f"#!/bin/bash\n#SBATCH -p compute\n#SBATCH -t {duration}\n#SBATCH --output=slurm-{jobName}.out\n#SBATCH --job-name={jobName}\n#SBATCH --mem={memory}G\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={cores}\n\n{codeLines}\n")

	return slurmFilename

####################################
# MAIN.
####################################

# Script parameters.
ucscName_human = "hg38" # Human is used as a point of reference in all comparisons.
minPCSsize     = 5 # number of base pairs. PCSs below this size will be ignored.

# Directory organization.
dirCode     = "/home/p/priscila-do/code" 
dirSlurm    = "/home/p/priscila-do/code/slurm"
dirMetadata = "/bucket/MillerU/Priscila/paper/data/inputs" # Directory "data" is available for download at Zenodo.
dirChains   = "/bucket/MillerU/Priscila/ucsc/chains"
dirFasta    = "/bucket/MillerU/Priscila/fasta"
dirPCSs     = "/flash/MillerU/Priscila/paper-validation/pcs"

# Slurm Parameters.
memory	   = 20  # GB
cores      = 1   # number of cores
duration   = f"2:00:00" # X:00:00 = X hours

print("******************************************************")
print("*   PCS (Perfectly conserved sequences) - CLUSTER    *")
print("* Extract PCSs from chains (ordered aligned blocks). *")
print("******************************************************")
print(f"- Directory for Python code: {dirCode}")
print(f"- Directory for Slurm scripts: {dirSlurm}")
print(f"- Directory for chain files (input): {dirChains}") # e.g.: *.hg38.all.chain.gz
print(f"- Directory for FASTA files (input): {dirFasta}")
print(f"- Directory for PCSs (output): {dirPCSs}")

# Check if all directories exist.
for dirPath in [dirCode, dirSlurm, dirMetadata, dirChains, dirFasta, dirPCSs]:
	if (not os.path.exists(dirPath)):
		print(f"ERROR! Directory not found ({dirPath}).")
		sys.exit()		

# Load species information.
speciesMetadataFilename = os.path.join(dirMetadata, "speciesMetadata.pickleProtocol4.pickle")
if(not os.path.isfile(speciesMetadataFilename)):
	print(f"ERROR! File with metadata for species not found ({speciesMetadataFilename}).")
	sys.exit()
speciesUCSCnames, divergenceTimes, commonNames = pickle.load(open(speciesMetadataFilename, 'rb'))
speciesUCSCnames = sorted(speciesUCSCnames, key=lambda species: (divergenceTimes[species], commonNames[species]))

# For each species.
for ucscName_other in speciesUCSCnames:

	# Check if input files exist.
	chainFilename = os.path.join(dirChains, f"{ucscName_other}.{ucscName_human}.all.chain.gz")
	dirFasta_qry  = os.path.join(dirFasta, ucscName_human, "chr")
	dirFasta_tar  = os.path.join(dirFasta, ucscName_other, "chr")

	if not(all(os.path.exists(filePath) for filePath in [chainFilename, dirFasta_qry, dirFasta_tar])):
		print(f"ERROR! {ucscName_other} : Check if chain file and FASTA files exist.\n\t- Chain filename: {chainFilename}\n\t- Fasta directory: {dirFasta_tar}\n\t- Fasta directory (human): {dirFasta_qry}")
		continue

	# Create slurm script.
	print(f"[{ucscName_other}] Create slurm script...")
	slurmFilename = os.path.join(dirSlurm,f"{ucscName_other}.pcs.slurm")
	createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, dirFasta_qry, dirFasta_tar, chainFilename, dirPCSs,
						ucscName_human, ucscName_other, minPCSsize)

	# Run slurm script.
	print(f"[{ucscName_other}] Running slurm script... {slurmFilename}")
	process = subprocess.Popen(['sbatch', slurmFilename])  # Replace with your command
	
