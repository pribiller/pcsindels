from collections import namedtuple
import time
import argparse
import os

"""  Data structure used to parse file from UCSC. 
""" 
Chain = namedtuple('Chain', 'score chainId tChrom tSize tStrand tStart tEnd qChrom qSize qStrand qStart qEnd')
Block = namedtuple('Block', 'igs tGap qGap')
Pcs   = namedtuple('Pcs', 'size tChrom tStrand tPosBeg qChrom qStrand qPosBeg')

"""  Functions to check input arguments. 
""" 
def dir_path(path):
	if os.path.isdir(path):
		return path
	else:
		raise argparse.ArgumentTypeError(f"ERROR! {path} is not a valid path. Check if the directory exists.")

def file_path(path):
	if os.path.isfile(path):
		return path
	else:
		raise argparse.ArgumentTypeError(f"ERROR! {path} is not a valid path. Check if the file exists.")


"""  Data structure to record usage of computational resources.
""" 
CompRes   = namedtuple('CompRes', 'time mem disk')

"""  Data structure to record computation time of scripts.
""" 
class Time:
	"""This class records the time taken for each step in the execution of a script.
	"""
	def __init__(self):
		self.times_lst = []

	def start(self):
		self.start_global = time.time()

	def stop(self):
		stop_global   = time.time()
		self.t_global = stop_global-self.start_global
		self.times_lst.append(("**Total time**",self.t_global))

	def startStep(self,desc,printTxt=True):
		self.step = desc
		self.start_local = time.time()
		if(printTxt): print(f"{self.step}...")

	def stopStep(self,printTxt=True):
		stop_local   = time.time()
		self.t_local = stop_local-self.start_local
		self.times_lst.append((self.step,self.t_local))
		if(printTxt): print(f"{self.step} : t={self.t_local}")

	def print(self):
		print("Stats on time:\n")

		t_str_lst = [f"{t:.2f}" for (desc,t) in self.times_lst]
		t_nbchars = max([len(t_str) for t_str in t_str_lst])+1
		d_nbchars = max([len(desc) for (desc,t) in self.times_lst])+1

		# Print title.
		print(f"{'='*d_nbchars}  {'='*t_nbchars}\n{'Step'.ljust(d_nbchars)}  {'Time (s)'.rjust(t_nbchars)}\n{'='*d_nbchars}  {'='*t_nbchars}")
		for (desc,t), t_str in zip(self.times_lst,t_str_lst):
			print(f"{desc.ljust(d_nbchars)}  {t_str.rjust(t_nbchars)}")
		print(f"{'='*d_nbchars}  {'='*t_nbchars}\n")
