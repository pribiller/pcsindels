"""Indels model. 
Implementation of the analytical equation described in the paper.
The function computes the expected PCS size distribution 
for each evolutionary time t.
"""

import numpy as np
from mpmath import * # zeta
import time
from kneed import DataGenerator, KneeLocator

mp.dps    = 30
mp.pretty = True

class IndelsSolver:

	# Constructor method to initialize the object
	def __init__(self, N, alpha, min_k):
		""" This class implements the indels model described in the paper,
			which is founded on a fragmentation system. In this evolutionary process,
			a long initial PCS fragments into smaller PCSs as time progresses, 
			causing mass loss and making base pairs within the PCSs 
			progressively more uncommon.
			
			:param N: Initial PCS size (in base pairs). It corresponds to the window size.
			:type desc: int
			:param alpha: 	It determines how frequent longer indels can occur. 
							The parameter alpha can take any value above 1: (1, ∞).
							If alpha is near 1, larger indels are more likely to occur.
							If alpha is above 5 (a hard upper limit set internally with [max_alpha]), 
							the model turns into a substitution only model.
			:type desc: float

			:param min_k: 	Minimum size of a PCS to be considered by the model (in base pairs). 
			:type desc: int

		"""

		# Model parameters.
		self.N	 = N
		self.alpha = alpha
		self.min_k = min_k

		# Pre-computed values.
		mu_k_alpha_all, Z_alpha = self.precomp_c_kt()
		self.mu_k_alpha_all = mu_k_alpha_all
		self.Z_alpha = Z_alpha

		self.sizes   = np.array(range(min_k,N+1))

	####################################
	# Indels model (analytical).
	def precomp_c_kt(self):
		"""This function pre-computes some terms of the analytical 
		equation that are invariant (regardless of the given evolutionary
		time t their value does not change). 
		"""
		max_alpha	  = 5
		is_subst_only = (self.alpha >= max_alpha)
		if(self.alpha in [1,2]):
			print(f"ERROR! Riemann zeta function undefined for α=1 (α and α-1 must be different from 1).")
			return (None, None)
		elif(is_subst_only):
			print(f"WARNING! High alpha (alpha={self.alpha}) : substitution only model will be used.")
			
		# Riemann zeta function for α and α-1.
		Z_alpha			= zeta(self.alpha)   if (not is_subst_only) else 1.0
		Z_alpha_minus_1 = zeta(self.alpha-1) if (not is_subst_only) else 1.0
		#print(f" Z(α)={Z_alpha}\n Z(α-1)={Z_alpha_minus_1}")
		def mu(k,alpha):
			return 2*((k-1)*Z_alpha + Z_alpha_minus_1)
			
		mu_k_alpha_all = [mu(k,self.alpha) for k in range(self.min_k,self.N+1)]
		return mu_k_alpha_all, Z_alpha
				
	def comp_c_kt(self,t):
		"""This function computes the PCS size distribution
		given an evolutionary time t.
		"""
		# Compute c(k,t).
		mu_N_alpha = self.mu_k_alpha_all[-1]
		c_kt	   = np.zeros(self.N-self.min_k+1)
		term_1	   = 4.0*t*self.Z_alpha
		term_2	   = 2.0*(t**2.0)*self.Z_alpha
		# 0<k<N
		for col, k in enumerate(range(self.min_k,self.N+1)):
			mu_k_alpha = self.mu_k_alpha_all[col]
			c_kt[col]  = mp.exp(-mu_k_alpha*t)*(term_1+term_2*(mu_N_alpha-mu_k_alpha))
		# k=N
		c_kt[-1] = mp.exp(-mu_N_alpha*t)*(1+4.0*t*self.Z_alpha)
		return c_kt
	
	def comp_c_kt_ln(self,t):
		# Compute c(k,t).
		mu_N_alpha = self.mu_k_alpha_all[-1]
		c_kt_ln	   = np.zeros(self.N-self.min_k+1)
		term_1	   = 4.0*t*self.Z_alpha
		term_2	   = 2.0*(t**2.0)*self.Z_alpha
		# 0<k<N
		for col, k in enumerate(range(self.min_k,self.N+1)):
			mu_k_alpha = self.mu_k_alpha_all[col]
			c_kt_ln[col]  = -mu_k_alpha*t + mp.log(term_1+term_2*(mu_N_alpha-mu_k_alpha))
		# k=N
		c_kt_ln[-1] = -mu_N_alpha*t + mp.log(1+4.0*t*self.Z_alpha)
		return c_kt_ln

	def find_ts(self,t_lims,max_nb_points):
		"""This function selects evolutionary times within an interval
		that are *sufficiently different* from one other, meaning that 
		they have distinct PCS size distributions. 
		
		These selected evolutionary times can be used later in a brute 
		force method, for example.
		
		*Sufficiently different* is defined based on two attributes:

			1. the sum of base pairs belonging to a PCS (which gives an idea about the mass loss); and
			2. the "knee" of the exponential curve (a point where the curve visibly bends, specifically from high slope to low slope).

		Run time is around 20 seconds.
		"""

		def comp_measure(t):
			# Warning! Other low ts also have elbow = min PCS size.
			if(t==0): 
				return (self.N, self.N)
			else: 
				c_kt     = self.comp_c_kt(t)
				totalBps = int(np.sum(self.sizes*c_kt))
				elbow	 = KneeLocator(self.sizes, c_kt, S=1.0, curve="convex", direction="decreasing").elbow if (totalBps > 0) else self.N*10
				return (elbow, totalBps)
					
		t_min, t_max  = t_lims
		queue   = [((t_min,comp_measure(t_min)),(t_max,comp_measure(t_max)))]
		all_ts  = []
		hist_ts = [t_min,t_max]
		all_val = []
		max_diff_val_addseltau = (5, 10) # Difference in elbow = 5 bps / Difference in total bps = 10 bps

		def add_t(info_t):
			t,val		    = info_t
			elbow, totalBps = val
			addVal = False
			if ((totalBps > 0) and (t not in all_ts)):
				if (len(all_val) < 1):
					addVal = True
				else:
					v_sel = min(all_val, key=lambda x: (abs(x[0]-elbow),abs(x[1]-totalBps)))
					if ((abs(v_sel[0]-elbow) >= max_diff_val_addseltau[0]) or (abs(v_sel[1]-totalBps) >= max_diff_val_addseltau[1])):
						addVal = True
			if (addVal):
				all_ts.append(t)
				all_val.append(val)
				
		def get_diff(info_t1,info_t2):
			return (abs(info_t1[1][0]-info_t2[1][0]),abs(info_t1[1][1]-info_t2[1][1]))
		
		def shouldAddQueue(info_t1,info_t2):
			diff_12=get_diff(info_t1,info_t2)
			return ((diff_12[0] >= max_diff_val_addseltau[0]) or (diff_12[1] >= max_diff_val_addseltau[1]))

		beg_time  = time.time()
		isTimeout = False
		timeout   = 300 # Max time allowed is 5 minutes. After that, it returns what got so far. Before: If there is nothing new in the last 6 minutes, abort. 
		while ((len(queue)>0) and (len(all_ts)<max_nb_points) and (not isTimeout)):
			# Get t interval.
			info_lb,info_ub  = queue.pop()
			t_lb = info_lb[0]
			t_ub = info_ub[0]
			
			# Compute new t.
			t_new = t_lb + (t_ub-t_lb)/2
			if (t_new in hist_ts):
				continue
			else:
				hist_ts.append(t_new)
			info_new = (t_new,comp_measure(t_new))
			
			# Add times only if there is not an added time with a similar PCS size distribution.
			# Similarity is measured based on the elbow and total bps values.
			add_t(info_lb)
			add_t(info_ub)

			# Queue only if the new time is "sufficiently" different from previous times.
			cur_time = time.time()
			if ((cur_time-beg_time) > timeout):
				isTimeout = True
				break

			if (shouldAddQueue(info_ub,info_new)): 
				queue.append((info_new, info_ub))
			if (shouldAddQueue(info_lb,info_new)): 
				queue.append((info_lb, info_new))
			# Sort time intervals in the queue giving priority to 
			# the intervals with the most different distribution.
			queue = sorted(queue, key=lambda x: get_diff(x[0],x[1]), reverse=True)
		end_time = time.time()

		extra_ts_nb = 100
		min_t = min([t for t in all_ts if t > 0])
		extra_ts = np.linspace(0, min_t, extra_ts_nb)[1:-1] # Exclude first and last points.

		all_ts   = sorted(set(list(extra_ts) + all_ts))
		# Check if there was a timeout.
		cur_time  = time.time()
		isTimeOut = ((cur_time-beg_time) > timeout)
		print(f"Select taus = {end_time-beg_time:.6f} seconds. (Timeout? {isTimeOut})")
		return all_ts
