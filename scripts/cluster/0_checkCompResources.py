import os
import re
import sys
import subprocess
import argparse
import time
import math
import itertools

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.dataset import Dataset
from utils.basicTypes import dir_path

def getComputationalResources(logFilename):
	"""This function extracts the time and the memory usage 
	from the information saved in the log.
	"""
	jobinival = 999999999
	jobid, time, memuse, disk_out, disk_tmp = jobinival, 0, 0, 0, 0

	with open(logFilename) as f:
		for line in f:
			m = re.match(r'Job ID: (\d+)', line)
			if m: jobid = m.group(1)
			m = re.match(r'Total disk usage: ([\d\.]+) GB', line)
			if m: disk_out = float(m.group(1))
			m = re.match(r'Temporary disk usage: ([\d\.]+) GB', line)
			if m: disk_tmp = float(m.group(1))
			m = re.match(r'Total time to run: ([\d\.]+) seconds', line)
			if m: time = float(m.group(1))
			m = re.match(fr"{jobid}.*\s+([\d\.]+)K", line)
			if m: memuse += float(m.group(1))

	# If memory was not saved, try to run sacct.
	if(memuse == 0):
		if (jobid != jobinival):
			# Define the command as a list of arguments
			command = ["sacct", "-j", jobid, "--format=CPUTime,MaxRSS"]
			# Run the command and capture the output.
			result = subprocess.run(command, capture_output=True, text=True)
			# Check if the command was successful.
			if (result.returncode == 0):
				for line in result.stdout.split("\n"):
					m = re.match(r".*\s+(\d+)K", line)
					if m: memuse += float(m.group(1))
			else:
				print(f"ERROR! {jobid=}, {time=}, {memuse=} in {logFilename}\n{result.stderr}")
		else:
			print(f"ERROR! {jobid=}, {time=}, {memuse=} in {logFilename}")
	return (jobid, time, memuse, disk_out, disk_tmp)

####################################
# MAIN.
####################################
# Usage:    python3 0_checkCompResources.py -script [script's name]
# Examples: python3 ~/code/cluster/0_checkCompResources.py -script extractPCS
#           python3 ~/code/cluster/0_checkCompResources.py -script computeWindows
#           python3 ~/code/cluster/0_checkCompResources.py -script setupEvolTimes
#           python3 ~/code/cluster/0_checkCompResources.py -script estimateEvolTimes
#           python3 ~/code/cluster/0_checkCompResources.py -script sampleEvolTimes
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Check usage of computational resources (memory and time).")
	parser.add_argument("-script",   help="Script to be checked.", type=str, required=True)

	args       = parser.parse_args()
	my_dataset = Dataset()
	
	scriptName = args.script
	print(f"Reading logs from script '{scriptName}'.")

	dictStr = ""
	inputs_rows = []
	inputs_cols = [None]
	if(scriptName == "extractPCS"):
		inputs_rows = my_dataset.speciesUCSCnames
	elif(scriptName == "computeWindows"):
		inputs_rows = my_dataset.chromLst
	elif(scriptName == "setupEvolTimes"):
		inputs_rows = [str(alpha) for alpha in my_dataset.alphas]
	elif(scriptName in ["estimateEvolTimes", "sampleEvolTimes"]):
		inputs_rows = my_dataset.speciesUCSCnames
		inputs_cols = [str(alpha) for alpha in my_dataset.alphas]

	name_nbchars    = max([len(str(inputStr)) for inputStr in inputs_rows])+1
	time_nbchars    = 9 # hh:mm:ss
	mem_nbchars     = 9
	disk_nbchars    = 8

	row_sepline_str = f"{'='*name_nbchars}"
	col_sepline_str = "  ".join([f"{'='*time_nbchars}  {'='*mem_nbchars}  {'='*disk_nbchars}" for i in inputs_cols])
	sepline_str     = f"{row_sepline_str}  {col_sepline_str}"
	row_title_str   = f"{'Desc.'.ljust(name_nbchars)}"
	col_title_str   = "  ".join([f"{'Time'.ljust(time_nbchars)}  {'Memory'.ljust(mem_nbchars)}  {'Disk'.ljust(disk_nbchars)}" for i in inputs_cols])
	title_str       = f"{row_title_str}  {col_title_str}"

	print(f"{sepline_str}\n{title_str}\n{sepline_str}")
	max_time = 86398     # 23:59:59
	max_mem  = 103809024 # 99 GB

	for inputStr_row in inputs_rows:
		rowStr = f"{str(inputStr_row).ljust(name_nbchars)}"
		for inputStr_col in inputs_cols:
			inputStr = (inputStr_row, inputStr_col) if (inputStr_col) else inputStr_row

			logFilename = ""
			if(scriptName == "extractPCS"):
				logFilename = my_dataset.getLogFilename_extractPCS(inputStr)
			elif(scriptName == "computeWindows"):
				logFilename = my_dataset.getLogFilename_computeWindows(inputStr)
			elif(scriptName == "setupEvolTimes"):
				logFilename = my_dataset.getLogFilename_setupEvolTimes(inputStr)
			elif(scriptName == "estimateEvolTimes"):
				logFilename = my_dataset.getLogFilename_estimateEvolTimes(*inputStr)
			elif(scriptName == "sampleEvolTimes"):
				logFilename = my_dataset.getLogFilename_sampleEvolTimes(*inputStr)

			if os.path.isfile(logFilename):
				jobid, t, memuse, diskGb, tmpGb = getComputationalResources(logFilename)
				
				# WARNING! It formats up to 24h (86399 seconds specifically).
				if (t >= 86399): print("WARNING! Job takes more than 24 hours!")

				######################################
				# Data info for documentation (exact values are used).

				# Adjust values.
				t      = t if (t > 0) else max_time
				memuse = int(memuse) if (memuse > 0) else max_mem
				diskGb = round(diskGb, 3) if (diskGb > 0) else 0
				tmpGb  = round(tmpGb, 3)  if (tmpGb > 0) else 0

				# Format values (hh:mm:ss and memory in GB)
				timeStr   = time.strftime('%H:%M:%S', time.gmtime(t)) if (t < 86399) else str(max_time)
				memuseGb  = math.ceil(memuse/1048576)
				diskGbStr = f"{diskGb:.3f}"

				rowStr += f"  {timeStr.ljust(time_nbchars)}  {(str(memuseGb)+'GB').ljust(mem_nbchars)}  {(diskGbStr+'GB').ljust(disk_nbchars)}"

				######################################
				# Data info for code. 
				# Ccomputational resource values are increased to create a bit of a slack in case 
				# the particular run was a lucky one (e.g., run on an exceptionally good computer).

				# Adjust values to add a bit of a slack (for code).
				t      = min(int(t*1.4),max_time) if (t > 0) else max_time
				memuse = min(int(memuse*1.15),max_mem) if (memuse > 0) else max_mem
				diskGb = round(diskGb*1.05, 3) if (diskGb > 0) else 0
				tmpGb  = round(tmpGb*1.05, 3) if (tmpGb > 0) else 0

				# Format values (hh:mm:ss and memory in GB)
				timeStr   = time.strftime('%H:%M:%S', time.gmtime(t)) if (t < 86399) else str(max_time)
				memuseGb  = math.ceil(memuse/1048576)
				diskGbStr = f"{diskGb:.3f}"

				# Save adjusted values for dictionary form.
				if(inputStr_row not in dictStr):
					dictStr += f'"{inputStr_row}": CompRes("{timeStr}", {memuseGb}, {diskGb}),'
			else:
				print(f"WARNING! Log file not found for input {inputStr}. Check if file exists: {logFilename}")
		print(rowStr)
		
	print(f"{sepline_str}\n")
	print(f"Printing in dictionary form:\n{{{dictStr[:-1]}}}")
