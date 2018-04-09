from gaussian_integral import *
import numpy as np

class CriticalOverlap:
	'''This is a supper class of gaussian integrals'''

	def __init__(self,LR,TF):	
		self.GI=GaussianIntegrals(LR,TF)
		#mean field equations
		self.overlap_eq = lambda del0,m: self.GI.Int_overlap_step_LR_qf(del0,m)
		self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_step_LR_qf(del0_c,m_c) 
		self.thres_m = 0.01

	def update_grids(self,dx,x_min,x_max):
		self.GI.update_grids(dx,x_min,x_max)
	
	def MFT_limit(self,method):
		if method==1: # einsum integrals  A,qf,qg,beta free
			self.overlap_eq = lambda del0,m: self.GI.Int_overlap_v2(del0,m)
			self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_v2(del0,m)
		elif method==2: # built-in integral A,qf,qg,beta free
			self.overlap_eq = lambda del0,m: self.GI.Int_overlap_v3(del0,m)
			self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_v3(del0,m)
		elif method==3: # limit A,p free and qf = qg
			self.overlap_eq = lambda del0,m: self.GI.Int_overlap_step_LR(del0,m)
			self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_step_LR(del0,m)
		elif method==4: # limit A,p and qf free
			self.overlap_eq = lambda del0,m: self.GI.Int_overlap_step_LR_qf(del0,m)
			self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_step_LR_qf(del0,m)
		elif method==5: # limit A->infty and p free 
			self.overlap_eq = lambda del0,m: self.GI.Int_overlap_large_A(del0,m)
			self.Int_TF2 = lambda del0,m: self.GI.Int_TF2_large_A(del0,m)

	#--------------Overlap------------------------------------------- 		
	def overlap(self,del0,m_init):
		m=m_init
		error=1.
		for i in range(50000):
			sol = self.overlap_eq(del0,m)
			error=np.abs(m-sol)
			m=sol
			#print m,del0,error,i
			if error<1e-7:
				return m
			if m<self.thres_m:
				return 0.
		return m
	
	#calculating the critical point
	def critical_point(self):
		#zero capacity
		m=10.
		m=self.overlap(0.1,m)
		if m<self.thres_m:
			return 0.1,0
		#coarse 100 grid
		d_del0=100.
		del0_0=(0.1,10000)
		del0_0,m=self.get_del0_init(m,del0_0,d_del0)
		
		if m<self.thres_m:
			return 0.1,0

		#coarse 10 grid
		d_del0=10.
		del0_0,m=self.get_del0_init(m,del0_0,d_del0)
		
		#coarse 1 grid
		d_del0=1.
		del0_0,m=self.get_del0_init(m,del0_0,d_del0)
		
		#fine grid
		d_del0=1e-1
		del0=np.arange(del0_0[0],del0_0[1],d_del0)
		m_last=m
		for l in range(len(del0)):
			m=self.overlap(del0[l],m)
			if m<self.thres_m: #	
				return del0[l-1],m_last
			else:
				m_last=m
		return del0[l],m
	
	#calculating the critical point
	def critical_point_grid(self):
		#zero capacity
		del0_min = 1e-3
		m=10.
		m=self.overlap(del0_min,m)
		#fine grid
		del0 = np.logspace(np.log10(del0_min),1,2000)
		m_last=m
		for l in range(1,len(del0)):
			m=self.overlap(del0[l],m)
			if m<self.thres_m: #	
				return del0[l-1],m_last
			else:
				m_last=m
		return del0[l],m
	
	def capacity(self):
		del0_c,m_c = self.critical_point_grid()
		if m_c<self.thres_m:
			return 0.
		else:	# prefactor limit
			alpha_c = del0_c/self.Int_TF2(del0_c,m_c)
			return alpha_c
		
	def get_del0_init(self,m,del0_0,d_del0):
		#medium 10 coarse grid	
		del0 = np.logspace(np.log10(0.05),0,300)
		m_last = m
		for j in range(len(del0)):
			m=self.overlap(del0[j],m)
			#print m,del0[j]
			if m<self.thres_m: 
				return (del0[j-1],del0[j]),m_last
			else:
				m_last=m	
		print 'Problem'
		return (del0_0[0],del0_0[1]),0.
	
	# self consistent equations	
	def self_consistent(self,del0_init,m_init,alpha):
		m = m_init
		del0 = del0_init
		error=1.
		for i in range(10000000):
			sol_m = self.overlap_eq(del0,m)
			sol_del0 = alpha *  self.Int_TF2(del0,m)
			error = np.abs(m-sol_m)
			if error<1e-10:
				m = sol_m
				del0 = sol_del0
				return del0,m
			if m<self.thres_m:
				return del0,0
			m = sol_m
			del0 = sol_del0
		return del0,m
