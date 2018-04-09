import numpy as np
from gaussian_integral import *
from critical_del0_m import *


class MFTCurves:
	'''This class provides the curves from the mft'''

	def __init__(self,LR,TF):
		self.GI=GaussianIntegrals(LR,TF)
		self.critical = CriticalOverlap(LR,TF)
		
	
	
	#computing capacity static mft
	def capacity_static_MFT(self):
		del0_c,m_c = self.critical.critical_point()
		if del0_c==0.1 and m_c<0.01:
			return 0.
		else:
			alpha_c = del0_c/(self.GI.myLR.gamma * self.GI.Int_TF2_v3(del0_c,m_c))
			return alpha_c
	

