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
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Check usage of computational resources (memory and time).")
	parser.add_argument("-script",   help="Script to be checked.", type=str, required=True)

	args       = parser.parse_args()
	my_dataset = Dataset()
	
	scriptName = args.script
	print(f"Reading logs from script '{scriptName}'.")

	dictStr = ""
	inputs  = []
	if(scriptName == "extractPCS"):
		inputs = my_dataset.speciesUCSCnames
	elif(scriptName == "computeWindows"):
		inputs = my_dataset.chromLst
	elif(scriptName == "setupEvolTimes"):
		inputs = [str(alpha) for alpha in my_dataset.alphas]
	elif(scriptName == "estimateEvolTimes"):
		inputs = [(UCSCname, str(alpha)) for alpha in my_dataset.alphas for UCSCname in my_dataset.speciesUCSCnames]

	name_nbchars = max([len(str(inputStr)) for inputStr in inputs])+1
	time_nbchars = 9 # hh:mm:ss
	mem_nbchars  = 9
	disk_nbchars = 8
	sepline_str  = f"{'='*name_nbchars}  {'='*time_nbchars}  {'='*mem_nbchars}  {'='*disk_nbchars}"
	print(f"{sepline_str}\n{'Desc.'.ljust(name_nbchars)}  {'Time'.ljust(time_nbchars)}  {'Memory'.ljust(mem_nbchars)}  {'Disk'.ljust(disk_nbchars)}\n{sepline_str}")
	max_time = 86398     # 23:59:59
	max_mem  = 103809024 # 99 GB

	for inputStr in inputs:

		logFilename = ""
		if(scriptName == "extractPCS"):
			logFilename = my_dataset.getLogFilename_extractPCS(inputStr)
		elif(scriptName == "computeWindows"):
			logFilename = my_dataset.getLogFilename_computeWindows(inputStr)
		elif(scriptName == "setupEvolTimes"):
			logFilename = my_dataset.getLogFilename_setupEvolTimes(inputStr)
		elif(scriptName == "estimateEvolTimes"):
			logFilename = my_dataset.getLogFilename_estimateEvolTimes(*inputStr)


		if os.path.isfile(logFilename):
			jobid, t, memuse, diskGb, tmpGb = getComputationalResources(logFilename)
			
			# WARNING! It formats up to 24h (86399 seconds specifically).
			if (t >= 86399): print("WARNING! Job takes more than 24 hours!")

			# Adjust values to add a bit of a slack.
			t      = int(t+t*0.5) if (t > 0) else max_time
			memuse = int(memuse + memuse*0.1) if (memuse > 0) else max_mem
			diskGb = round(diskGb + diskGb*0.05, 2) if (diskGb > 0) else 0
			tmpGb  = round(tmpGb  + tmpGb*0.05, 2) if (tmpGb > 0) else 0

			# Format values (hh:mm:ss and memory in GB)
			timeStr   = time.strftime('%H:%M:%S', time.gmtime(t)) if (t < 86399) else str(max_time)
			memuseGb  = math.ceil(memuse/1048576)
			diskGbStr = f"{diskGb:.3f}"
			print(f"{str(inputStr).ljust(name_nbchars)}  {timeStr.ljust(time_nbchars)}  {(str(memuseGb)+'GB').ljust(mem_nbchars)}  {(diskGbStr+'GB').ljust(disk_nbchars)}")

			dictStr += f'"{inputStr}": CompRes("{timeStr}", {memuseGb}, {diskGb}),'
		else:
			print(f"WARNING! Log file not found for input {inputStr}. Check if file exists: {logFilename}")
	
	print(f"{sepline_str}\n")
	print(f"Printing in dictionary form:\n{{{dictStr[:-1]}}}")
