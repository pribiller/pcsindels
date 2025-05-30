"""Indels model. 
Implementation of the analytical equation described in the paper.
The function computes the expected PCS size distribution 
for each evolutionary time t.
"""

import numpy as np
from mpmath import * # zeta
import time

from collections import namedtuple, defaultdict

from kneed import DataGenerator, KneeLocator
from scipy.optimize import lsq_linear
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

import sys

mp.dps    = 30
mp.pretty = True

# Ways to summarize the PCS size distribution in few measures.
FuncExp = namedtuple('FuncExp', 'elbow totalBps')
FuncLin = namedtuple('FuncLin', 'slope intercept')

####################################################################
# Useful functions for estimation method.

# Linearize PCS size distribution.
def linearize(PCSsizes, PCScnts):
	y = np.log(PCScnts + [1e-1])   # "real" y-axis.
	x = PCSsizes+[max(PCSsizes)+1] # "real" x-axis.
	# Prepare the design matrix A (with a column of ones for the intercept)
	A = np.vstack([np.ones(len(x)), x]).T  # Shape (100, 2)
	# Target vector b
	b = y  # Shape (100,)
	# Define bounds for the coefficients
	# Ideally, the slope should be negative, but in low samples this might not be the case.
	bounds = ([-np.inf, -np.inf], [np.inf, -1e-10])  # Intercept can be any value, slope <= 0
	# Solve the least squares problem
	result = lsq_linear(A, b, bounds=bounds, tol=1e-20)
	# Extract the optimized parameters
	intercept_opt, slope_opt = result.x
	return (slope_opt,intercept_opt)

def getSampleSize(PCSsizeDistrib):
	return sum(list(PCSsizeDistrib.values()))

####################################################################
# Sampling related methods.
def samplePCSsNonPar(input):
	LocalProcRandGen = np.random.RandomState()
	PCSsizes, PCSsizes_probs, N = input
	sample = defaultdict(int)
	N_samp = 0
	while(N_samp < N):
		sampleSize = min(N-N_samp,10000)
		PCSsizes_sampled = LocalProcRandGen.choice(PCSsizes, sampleSize, p=PCSsizes_probs, replace=True)
		for PCSsize in PCSsizes_sampled: sample[PCSsize] += 1
		N_samp += sampleSize
	return sample
	
def samplePCSs(PCSsizes, PCSsizes_probs, N, nbCores=1):
	if(nbCores == 1):
		return samplePCSsNonPar((PCSsizes, PCSsizes_probs, N))
	else:
		# Prepare parallel inputs.
		N_par = int(N/(nbCores))
		N_rem = N % (nbCores)
		parallelInputs = [(PCSsizes, PCSsizes_probs, N_par) for idxSample in range(nbCores-N_rem)]
		parallelInputs.extend([(PCSsizes, PCSsizes_probs, N_par+1) for idxSample in range(N_rem)])		
		# Parallel sampling.
		sample = {}
		with futures.ProcessPoolExecutor(nbCores) as pool:
			for result in pool.map(samplePCSsNonPar, parallelInputs):
				PCSsizes_sampled, PCScnts_sampled = result
				# Aggregating results.
				for (PCSsize,PCScnt) in zip(PCSsizes_sampled,PCScnts_sampled):
					if (PCSsize not in sample.keys()): sample[PCSsize] = 0
					sample[PCSsize] += PCScnt
		return sample


###########################################################
# Evolutionary model with one parameter α to indicate the 
# propensity of insertions and deletions to occur.
class IndelsModel:

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


####################################################################
# Class to estimate evolutionary times.	
class EvolTimesSolver:
	
	def __init__(self, ts):
		self.ts     = ts
		self.ts_KDE = []

	def computePcsProbs_t(self, t, model):
		""" WARNING! Sample bias.
			The solver ignores that some evolutionary times might have considerably 
			more PCSs **under** ``minPCSsize``, thus creating a sample with PCSs 
			above ``minPCSsize`` would be very unlikely.
		"""
		# Computes the expected counts for each PCS size given an evolutionary time t.
		PCSprobs_t = model.comp_c_kt(t) # idx 0 = min PCS size
		PCSprobs_t[-1] = 0.0 # Ignore max PCS size (=window size).
		# Normalize counts to get probabilities.
		totalCount = sum(PCSprobs_t)
		PCSprobs_t = PCSprobs_t/totalCount
		return (model.sizes, PCSprobs_t)

	def computeSample_t(self, t, model, nbDistribs_sampled, nbPCSs_sampled, minDiffPcsSizes):
		# Consistency check.
		minDiffPcsSizes = max(minDiffPcsSizes, 1)
		if(nbPCSs_sampled < minDiffPcsSizes):
			print(f"ERROR! Number of PCSs sampled (={nbPCSs_sampled}) is smaller than the minimum number of PCSs needed (={minDiffPcsSizes}).")
			sys.exit()
		t_PCSsizes, t_PCSprobs = self.computePcsProbs_t(t, model)
		sample = []
		for sampleIdx in range(nbDistribs_sampled):
			# Sample PCS sizes from the PCS size distribution. 
			# Use same conditions as the ones found in observed data (same count of PCSs).
			# Consider only PCS sizes above min. PCS size.
			PCSdistrib_samp = {}
			nb_attempts	    = 0
			nb_attempts_max = 100
			# Discard PCS size distribution if number of distinct PCS sizes is below [minDiffPcsSizes].
			while (len(PCSdistrib_samp) < minDiffPcsSizes): 
				PCSdistrib_samp = samplePCSs(t_PCSsizes,t_PCSprobs, nbPCSs_sampled, 1)
				# Makes the sampling process easier if it is taking too much time.
				nb_attempts += 1
				if(nb_attempts > nb_attempts_max): 
					nb_attempts = 0
					minDiffPcsSizes=minDiffPcsSizes-1
			sample.append(PCSdistrib_samp)
		return sample

	def kdeInputs(self, PcsSizeDistrib_lst):
		""" This function handles the information from a list of
			windows to be used as input for the inference method. 
		"""
		# Exponential behavior is expected to be linear in (constrained) log-y.
		slopes     = []
		intercepts = []
		for PcsSizeDistrib in PcsSizeDistrib_lst:
			PCSsizes = list(sorted(PcsSizeDistrib.keys()))
			PCScnts  = [PcsSizeDistrib[PCSsize] for PCSsize in PCSsizes]
			slope, intercept = linearize(PCSsizes, PCScnts)
			slopes.append(slope)
			intercepts.append(intercept)

		# If variables are correlated, gaussian_kde gets stuck.
		# gaussian_kde does not currently support data that lies in a 
		# lower-dimensional subspace of the space in which it is expressed. 

		# Combine x and y into a 2D array
		data     = np.vstack([slopes, intercepts]) # Gaussian KDE: shape (# of dims, # of data).
		data_pca = data.T # PCA: shape (n_samples, n_features) (flipped shape compared to Gaussian KDE).
		return data, data_pca

	def computePcsKde_t(self, t_samps):
		data, data_pca = self.kdeInputs(t_samps)
		kde  = None
		pca  = None
		try:
			kde  = gaussian_kde(data, bw_method='scott') # 2-D array with shape (# of dims, # of data).
		# This exception was happening for the slope and nb. bps; they were not 2 distinct dimensions in few cases.
		# My guess is that the slope becomes random/uninformative when there are only few points.
		# In our case, this exception happens usually for the smallest sample size allowed (sampleSizeRef=5),
		# and rarely for other sample sizes (sampleSizeRef=8,11).
		except np.linalg.LinAlgError as e:
			print(f"[{t=} {nbSamplesPerTau=} {sampleSize=}] Gaussian KDE exception caught: {e}\n{slopes}\n{intercepts}")
			# Perform PCA to reduce dimensions
			pca = PCA(n_components=1)  # Reduce to 1 dimension
			pca.fit(data_pca)
			reduced_data = pca.transform(data_pca)
			# Now try Gaussian KDE on the reduced data.
			kde = gaussian_kde(reduced_data.T) # Gaussian KDE: shape (# of dims, # of data)
		pdf_values = kde.evaluate(data) if (pca == None) else kde.evaluate(pca.transform(data_pca).T)
		return (kde, pca, pdf_values.max())

	def initialize_ts(self, model, nbDistribs_sampled, nbPCSs_sampled, minDiffPcsSizes):
		self.ts_KDE = []
		for t in self.ts:
			samp_info = self.computeSample_t(t, model, nbDistribs_sampled, nbPCSs_sampled, minDiffPcsSizes)
			kde_t, pca_t, max_likelihood_t = self.computePcsKde_t(samp_info)
			self.ts_KDE.append((kde_t, pca_t, max_likelihood_t))

	def estimate_ts(self, pcsSizeDistrib_lst):
		""" This function estimates the evolutionary times of ``N``
			PCS size distributions specified in ``pcsSizeDistrib_lst``.

		"""
		# Preparing input data.
		data, data_pca = self.kdeInputs(pcsSizeDistrib_lst)
		# Compute the posterior for each window.
		nbrows	  = len(self.ts)
		nbcols	  = len(pcsSizeDistrib_lst)
		M_posterior = np.zeros((nbrows,nbcols))
		for t_idx, (t, t_KDE) in enumerate(zip(self.ts, self.ts_KDE)):
			kde_t, pca_t, max_likelihood_t = t_KDE
			# Compute likelihood that observed data was sampled from the same distribution.
			likelihood_obs_t   = kde_t.pdf(data) if (pca_t == None) else kde_t.pdf(pca_t.transform(data_pca).T)
			M_posterior[t_idx] = likelihood_obs_t/max_likelihood_t
		return M_posterior
	
####################################################################
# Class to select evolutionary times.
class EvolTimesSelector:

	# Declare the named tuple as a class attribute
	_EvolTime = namedtuple('EvolTime', ['t', 'desc'])

	def __init__(self, model):
		self.model = model
		self.ts    = []
		# Criteria to determine sufficiently different functions:
		# 1. Difference in elbow >= 5 bps
		# 2. Difference in total bps >= 10 bps
		self.criteriaDiffFuncs = FuncExp(5, 10)
		# Max time allowed to find evolutionary times is 5 minutes. 
		# After that, it returns what was found so far.
		self.find_ts_timeout = 300 # seconds
		self.find_ts_begTime = 0
		# Number of extra evolutionary times added between (0, min_t).
		self.nb_ts_extra = 100

	def isTimeout(self):
		return ((time.time()-self.find_ts_begTime) > self.find_ts_timeout)

	def get_ts(self):
		return sorted([t_info.t for t_info in self.ts])

	def new_t(self, t):
		"""This function instantiates a new evolutionary time object.
		This object has two attributes:

			1. ``t`` : The time itself;
			2. ``desc`` : A *summarized description* of the expected PCS size distribution at time ``t``.
		
		The PCS size distribution is characterized by two measures:
			1. The *elbow* of the PCS size distribution function;
			2. The total number of base pairs which are PCSs.

		"""
		N     = self.model.N
		sizes = self.model.sizes
		# Warning! Other low ts also have elbow = min PCS size.
		if(t==0): 
			return self._EvolTime(t, FuncExp(N, N))
		else: 
			c_kt     = self.model.comp_c_kt(t)
			totalBps = int(np.sum(sizes*c_kt))
			elbow	 = KneeLocator(sizes, c_kt, S=1.0, curve="convex", direction="decreasing").elbow if (totalBps > 0) else N*10
			return self._EvolTime(t, FuncExp(elbow, totalBps)) 

	def funcDiff(self,t_1,t_2):
		return FuncExp(abs(t_1.desc.elbow-t_2.desc.elbow), abs(t_1.desc.totalBps-t_2.desc.totalBps))

	def isSufficientlyDiff(self,info_t1,info_t2):
		diff=self.funcDiff(info_t1,info_t2)
		return ((diff.elbow >= self.criteriaDiffFuncs.elbow) or (diff.totalBps >= self.criteriaDiffFuncs.totalBps))

	def add_t(self, new_t):
		"""This function appends an evolutionary time to the list 
		if the new value *significantly differs* from the existing times.

		"""
		if (new_t.desc.totalBps > 0):
			# Check if all elements are True
			all_diff = all([self.isSufficientlyDiff(new_t, old_t) for old_t in self.ts])
			if(all_diff): self.ts.append(new_t)
			
	def add_ts(self,min_t,max_t,nb_ts):
		""" This function adds N points (if they are not already added)
			to the list of selected evolutionary times. These points are 
			linearly distributed in the interval ``(min_t, max_t)``.
		"""
		ts_all   = self.get_ts()
		extra_ts = np.linspace(min_t,max_t,nb_ts)[1:-1] # Exclude first and last points.
		new_ts   = set(extra_ts) - set(ts_all)
		self.ts.extend([self.new_t(t) for t in new_ts])

	def find_ts(self,t_lims,max_nb_points):
		"""This function selects evolutionary times within an interval
		that are *sufficiently different* from one other, meaning that 
		they have distinct PCS size distributions. 
		
		These selected evolutionary times can be used later in a brute 
		force method, for example.
		
		*Sufficiently different* is defined based on two attributes:

			1. the sum of base pairs belonging to a PCS (which gives an idea about the mass loss); and
			2. the "knee" of the exponential curve (a point where the curve visibly bends, specifically from high slope to low slope).

		Run time is around 25 seconds.
		"""
		
		# Initialize parameters.
		t_min, t_max  = t_lims
		queue   = [(self.new_t(t_min),self.new_t(t_max))]
		hist_ts = [t_min,t_max]

		# Controls timeout.
		self.find_ts_begTime = time.time()
		while ((len(queue)>0) and (len(self.ts)<max_nb_points) and (not self.isTimeout())):
			# Get t interval.
			info_lb,info_ub  = queue.pop()
			t_lb = info_lb.t
			t_ub = info_ub.t
			
			# Compute new t.
			t_new = t_lb + (t_ub-t_lb)/2
			if (t_new in hist_ts):
				continue
			else:
				hist_ts.append(t_new)
			info_new = self.new_t(t_new)
			
			# Add times only if there is not an added time with a similar PCS size distribution.
			# Similarity is measured based on the elbow and total bps values.
			self.add_t(info_lb)
			self.add_t(info_ub)

			# Update queue if there is some computation time left.
			if(not self.isTimeout()):
				# Queue only if the new time is "sufficiently" different from previous times.
				if (self.isSufficientlyDiff(info_ub,info_new)):  queue.append((info_new, info_ub))
				if (self.isSufficientlyDiff(info_lb, info_new)): queue.append((info_lb, info_new))
				# Sort time intervals in the queue giving priority to 
				# the intervals with the most different distribution.
				queue = sorted(queue, key=lambda t_interv: self.funcDiff(t_interv[0], t_interv[1]), reverse=True)

		# Add a couple of extra evolutionary times.
		ts_so_far   = self.get_ts()
		self.add_ts(0, min([t for t in ts_so_far if t > 0]), self.nb_ts_extra)
		return self.get_ts()
