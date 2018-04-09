from gaussian_integral import *
import numpy as np

class CriticalOverlap:
	'''This is a supper class of gaussian integrals'''

	def __init__(self,LR,TF):	
		self.gaussian_integral=GaussianIntegrals(LR,TF)
		
	def update_grids(self,dx,x_min,x_max):
		self.gaussian_integral.update_grids(dx,x_min,x_max)

	#--------------Overlap------------------------------------------- 		
	def overlap(self,del0,m_init):
		m=m_init
		error=1.
		for i in range(1000000):
			sol=self.gaussian_integral.Int_overlap(del0,m)
			error=np.abs(m-sol)
			m=sol
			if error<1e-7:
				return m
			if m<1e-2:
				return 0.
		return m
