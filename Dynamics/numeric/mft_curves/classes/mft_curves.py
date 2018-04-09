import numpy as np
from scipy.optimize import brentq
from gaussian_integral import *
from critical_del0_m import *


class MFTCurves:
	'''This class provides the curves from the mft'''

	def __init__(self,LR,TF):
		self.GI=GaussianIntegrals(LR,TF)
		self.critical = CriticalOverlap(LR,TF)
		
		# if del0_c or m_c is calculated
		self.del0_c = 0.
		self.m_c = 0.
		self.del1 = 0. # int critical del1
	

	#overlap curve static mft
	def overlap_static_MFT(self,alpha,q_init=10.,return_m=True):
		# defining the model
		myq = lambda y:self.critical.overlap(y,q_init)
		myfield = lambda x:alpha*self.GI.myLR.gamma*self.GI.Int_TF2(x,myq(x))-x
		del0 = brentq(myfield,0,2000.)
		q = myq(del0)
		mean_TF = self.GI.Int_TF(del0,q)
		mean_TF2 = self.GI.Int_TF2(del0,q)
		var_TF = mean_TF2 - mean_TF**2
		m = q/np.sqrt(self.GI.myLR.intg2 * var_TF)
		print 'alpha=',alpha,'overlap=',m,'del0=',del0
		return np.array([(m,q),del0])

	

