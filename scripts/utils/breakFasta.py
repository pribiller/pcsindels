import os
import argparse
import re
import sys

def file_path(filepath):
	if os.path.isfile(filepath):
		return filepath
	else:
		raise argparse.ArgumentTypeError(f"readable_file: '{filepath}' is not a valid file.")

def dir_path(dirpath):
	if os.path.exists(dirpath):
		return dirpath
	else:
		raise argparse.ArgumentTypeError(f"readable_dir: '{dirpath}' is not a valid directory.")

####################################
# MAIN.
####################################
# Usage example: python3 breakFasta.py -fastaDir /media/pri/Pics/IGS/miniFasta -prefix test -size 20
# python3 breakFasta.py -fastaDir /path/to/fasta/mm39-Mouse-Jun2020 -prefix mm39 -size 20000

parser = argparse.ArgumentParser(description="Break FASTA")

parser.add_argument("-fastaDir", help="directory of FASTA file.", type=dir_path, required=True)
parser.add_argument("-prefix",   help="prefix of FASTA file.", required=True)
parser.add_argument("-size",     help="number of BPs (base pairs/nucleotides) per file", type=int, required=True)
args = parser.parse_args()

fastaDir    = args.fastaDir
fastaPrefix = args.prefix
fastaSize   = args.size
fastaFileName = os.path.join(fastaDir,fastaPrefix+".fa")

print("*****************************************************")
print("*                   Break FASTA                     *")
print("*****************************************************")

# Create directories.
print("Creating directories...")
with open(fastaFileName) as f:
	for line in f:
		#>chr1
		m = re.match(r">([\S]+)\s*$", line)
		if m:
			chrName = m.group(1)
			chrDir  = os.path.join(fastaDir,"chr",chrName)
			if not os.path.exists(chrDir):
				# Create a new directory because it does not exist
				os.makedirs(chrDir)

# Break FASTA into directories.
print("Reading FASTA files...")
with open(fastaFileName) as f:
	chrName = ""
	chrDir  = ""
	dna     = ""
	absPos  = 0
	begPos  = -1
	for line in f:
		#>chr1
		m = re.match(r">([\S]+)\s*$", line)
		if m:
			
			# Save last fragment from previous chromosome.
			if (len(dna) > 0):
				#print("absPos=" + str(begPos) + " / end=" + str(begPos+len(dna)))
				chrFilename = os.path.join(chrDir, fastaPrefix + "." + chrName + ".b" + str(begPos) + "e" + str(begPos+len(dna)) + ".fa")
				with open(chrFilename, "w") as dna_file:
					dna_file.write(dna)

			# Start a new chromosome.
			chrName = m.group(1)
			chrDir  = os.path.join(fastaDir,"chr",chrName)
			print(chrName)
			dna     = ""
			absPos  = 0
			begPos  = 0
		else:
			m = re.match(r"([\S]+)\s*$", line)
			if m:
				dnaFragment = m.group(1)
				dnaFragSize = len(dnaFragment)
				curSize     = len(dna)

				while ((curSize+dnaFragSize) >= fastaSize):
					# Save fragment.
					dna    += dnaFragment[:(fastaSize-curSize)]
					absPos += len(dnaFragment[:(fastaSize-curSize)])
					# genomeName.chrName.bXXXeXXX.fa
					#print("absPos=" + str(begPos) + " / end=" + str(begPos+fastaSize) + " / fastaSize=" + str(fastaSize))
					chrFilename = os.path.join(chrDir, fastaPrefix + "." + chrName + ".b" + str(begPos) + "e" + str(begPos+fastaSize) + ".fa")
					with open(chrFilename, "w") as dna_file:
						dna_file.write(dna)

					# Update DNA fragment.
					dnaFragment = dnaFragment[(fastaSize-curSize):]
					dnaFragSize = len(dnaFragment)

					# Start new fragment.
					dna     = ""
					begPos  = absPos
					curSize = 0

				dna    += dnaFragment
				absPos += dnaFragSize

	# Save last fragment.
	if (len(dna) > 0):
		#print("absPos=" + str(begPos) + " / end=" + str(begPos+len(dna)))
		chrFilename = os.path.join(chrDir, fastaPrefix + "." + chrName + ".b" + str(begPos) + "e" + str(begPos+len(dna)) + ".fa")
		with open(chrFilename, "w") as dna_file:
			dna_file.write(dna)
