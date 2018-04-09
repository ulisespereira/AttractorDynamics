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
	#calculating the critical point
	def critical_point(self,del0_init):
		d_del0=1e-2
		del0=np.arange(del0_init,10000.,d_del0)
		m=100.
		m_last=m
		for i in range(len(del0)):
			m=self.overlap(del0[i],m)
			if m<0.01 and 0<i: #
				return (del0[i-1],m_last)
			elif m<0.01 and i==0:
				return (0,0)
			else:
				m_last=m
			print del0[i],m
		return del0[i],m
