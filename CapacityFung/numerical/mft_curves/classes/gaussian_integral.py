from transfer_function import *
from learning_rule import *
import numpy as np
from scipy.stats import multivariate_normal as mvnormal
import time
class GaussianIntegrals:
	'''This is a supper class of gaussian integrals'''

	def __init__(self,LR,TF):
		
		self.myTF=TF # a transfer function class
		self.myLR=LR # learning rule class
		
		# integration variables
		self.dx=0.1 #
		x_min=-10.
		x_max=10.
		self.xgrid = np.arange(x_min,x_max,self.dx) #grid on x
		self.normal_pdf =self.myLR.std_normal(self.xgrid) #standard normal pdf
		self.eta=self.myTF.TF(self.xgrid)# trasnformation tf
		self.f_eta=self.myLR.f(self.eta)# trasnformtion eta by f
		self.g_eta=self.myLR.g(self.eta)# trasnformtion eta by g
		# defining used vectors/matrices

	def update_grids(self,dx,x_min,x_max):
		# integration variables
		self.dx=dx #
		self.xgrid = np.arange(x_min,x_max,self.dx) #grid on x
		self.normal_pdf =self.myLR.std_normal(self.xgrid) #standard normal pdf
		self.eta=self.myTF.TF(self.xgrid)# trasnformation tf
		self.f_eta=self.myLR.f(self.eta)# trasnformtion eta by f

	# integral corresponding to the overlap	
	def Int_overlap(self,del0,m):	
		meanfield=self.myLR.Amp*self.f_eta*m
		sdfield=np.sqrt(del0)*self.xgrid
		steady_state=self.myTF.TF(np.add.outer(meanfield,sdfield))
		g=self.normal_pdf*self.g_eta
		sol=self.dx*self.dx*np.einsum('i,ij,j',g,steady_state,self.normal_pdf)
		return sol
	
	# integral TF del1=del0 static MFT
	def Int_TF(self,del0,m):
		meanfield = self.myLR.Amp*self.f_eta*m
		sdfield = np.sqrt(del0)*self.xgrid
		tf = self.myTF.TF(np.add.outer(meanfield,sdfield))
		sol=self.dx*self.dx*np.einsum('i,ij,j',self.normal_pdf,tf,self.normal_pdf)
		return sol

	# integral TF**2 del1=del0 static MFT
	def Int_TF2(self,del0,m):
		meanfield = self.myLR.Amp*self.f_eta*m
		sdfield = np.sqrt(del0)*self.xgrid
		tf = self.myTF.TF(np.add.outer(meanfield,sdfield))
		sol=self.dx*self.dx*np.einsum('i,ij,ij,j',self.normal_pdf,tf,tf,self.normal_pdf)
		return sol


