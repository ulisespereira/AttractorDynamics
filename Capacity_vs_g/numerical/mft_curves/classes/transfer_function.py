import numpy as np
from scipy.optimize import brentq


class TransferFunction:
	'''this class is for different transfer functions'''

	def __init__(self,paramTF):

		self.which_tf=paramTF[0]
		self.rm=paramTF[1]
		self.b=paramTF[2]
		self.h0=paramTF[3]	
		if self.which_tf=='tanh':
			self.q_tanh=paramTF[4]

	
	#TF
	def TF(self,h):
		# sigmoidal TF
		if self.which_tf=='sig':
			phi = self.rm/(1.+np.exp(-self.b * (h-self.h0)))
		# tanh TF
		elif self.which_tf=='tanh':
			phi =  (self.rm/2.) * (2 * self.q_tanh-1. + np.tanh(self.b * (h - self.h0)))
		return phi
	#DTF
	def DTF(self,h):
		#sigmoidal
		if self.which_tf=='sig':
			phi=(self.rm * np.exp(-self.b*(h-self.h0)) * self.b) / ((1.+np.exp(-self.b*(h-self.h0)))**2)
		#tanh
		elif self.which_tf=='tanh':
			phi = 0.5 * self.b * self.rm * (1./np.cosh(self.b * (h - self.h0))**2)
		return phi
	#D2TF
	def D2TF(self,h):
		if self.which_tf=='sig':
			phi1 = 2*(self.b/self.rm)*self.TF(h)*self.DTF(h)*np.exp(-self.b*(h-self.h0))
			phi2 = -((self.b**2)/self.rm)*(self.TF(h)**2)*np.exp(-self.b*(h-self.h0))
			phi = phi1 + phi2
		elif self.which_tf=='tanh':
			phi = -(self.b**2) * self.rm * (1./np.cosh(self.b * (h - self.h0))**2) * np.tanh(self.b * (h - self.h0))
		return phi
        
	# Int TF
	def IntTF(self,h):
		if self.which_tf=='sig':
			#phi = self.rm*(h + (1/self.b)*np.log( (1+np.exp(self.b*(self.h0-h)))/(1+np.exp(self.b*self.h0)) ) )# Yonatan's choice	
			#phi=self.rm*(h+np.log(1+np.exp(-self.b*(h-self.h0)))/self.b-self.h0)
			phi = self.rm * (h + (1/self.b)*np.log(1 + np.exp(-self.b * (h -self.h0))))
		elif self.which_tf=='tanh':
			phi = (self.rm/2)*( (2 * self.q_tanh -1.) * h + np.log(np.cosh(self.b * (h - self.h0)))/self.b)
		return phi
	
	#inverse TF
	def TF_Inv(self,r):
		if self.which_tf=='sig':
			phi=self.h0 + (1./self.b) * np.log(r/(self.rm-r))
		elif self.which_tf=='tanh':
			phi = self.h0 + (1./self.b) * np.arctanh(2*r/self.rm - (2 * self.q_tanh -1))
		return phi
	
	#derivative inverse TF
	def DTF_Inv(self,r):
		if self.which_tf=='sig':
			phi = (1./self.b) * (self.rm/(r * (self.rm-r)))
		elif self.which_tf=='tanh':
			y=2*r/self.rm - (2 * self.q_tanh -1)
			phi = 2./(self.b * self.rm * (1 - y**2))
		return phi
