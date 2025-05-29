"""

Output states
^^^^^^^^^^^^^
==========   ============  =================
UCSC name        Size        Nb. of files 
==========   ============  =================
panPan3           2.84 GB             154815
panTro6           2.84 GB             154808
gorGor6           2.84 GB             155144
ponAbe3           2.85 GB             156067
papAnu4           2.76 GB             203591
macFas5           2.74 GB             153166
rhiRox1           2.70 GB             273572
chlSab2           2.60 GB             140996
nasLar1           2.81 GB             453159
rheMac10          2.77 GB             150086
calJac4           2.70 GB             145394
tarSyr2           3.22 GB             478822
micMur2           2.27 GB             129366
galVar1           2.97 GB             310790
mm39              2.54 GB             136442
oryCun2           2.55 GB             138594
rn7               2.47 GB             132485
vicPac2           2.02 GB             378096
bisBis1           2.75 GB             577597
felCat9           2.35 GB             128505
manPen1           2.05 GB             185418
bosTau9           2.53 GB             136904
canFam6           2.15 GB             115734
musFur1           2.25 GB             127030
neoSch1           2.24 GB             126686
equCab3           2.33 GB             128329
myoLuc2           1.89 GB             108148
susScr11          2.33 GB             125395
enhLutNer1        2.26 GB             167835
triMan1           2.89 GB             160112
macEug2           2.86 GB             367992
ornAna2           1.86 GB             272083
aptMan1           1.42 GB              95522
galGal6           0.99 GB              53546
thaSir1           1.33 GB              76505
aquChr2           1.11 GB              60432
melGal5           1.05 GB             280828
xenLae2           2.53 GB             241082
xenTro10          1.35 GB              72683
danRer11          1.56 GB              84923
Total        **92.57 GB**  **7538682 files**
==========   ============  =================

"""
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
