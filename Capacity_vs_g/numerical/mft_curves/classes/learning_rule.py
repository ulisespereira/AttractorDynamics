import numpy as np
from transfer_function import *
from scipy import integrate
from scipy.optimize import brentq
	
	
	
	

class LearningRule:
	''' This class gives the learning rule'''
	def __init__(self,paramsLR,TF):

		self.myTF=TF
		#parameters function g and f
		self.xf=paramsLR[0]
		self.xg=paramsLR[1]
		self.betaf=paramsLR[2]
		self.betag=paramsLR[3]
		self.qf=paramsLR[4]
		self.Amp=paramsLR[5]
		
		# here it is very important do it in this order	
		self.qg=self.Qg()
		self.intg2=self.Eg2()
		self.intf2=self.Ef2()	
		#print 'Amp',self.Amp,myint
		self.gamma=self.intf2*self.intg2*self.Amp*self.Amp
		
	
	def std_normal(self,x):
		sigma=1. # standard normal random variable passed thorugh transfer functiion
		mu=0
		pdf=(1./np.sqrt(2 * np.pi * sigma**2))*np.exp(-(1./2.)*((x-mu)/sigma)**2)
		return pdf
	
	# separable functions learning process leanring rule f and g
	def f(self,x):
		return 0.5*(2*self.qf-1.+np.tanh(self.betaf*(x-self.xf)))

	def g(self,x):
		return 0.5*(2*self.qg-1.+np.tanh(self.betag*(x-self.xg)))

	def Qg(self):# mean of f**2
		return brentq(self.Eg,0.,1.)

	def Eg(self,q):# mean of g
		self.qg=q
		fun=lambda x:self.std_normal(x)*self.g(self.myTF.TF(x))
		var,err=integrate.quad(fun,-10.,10.)
		return var

	def Ef2(self):# mean of f**2
		fun=lambda x:self.std_normal(x)*self.f(self.myTF.TF(x))*self.f(self.myTF.TF(x))
		var,err=integrate.quad(fun,-10.,10.)
		return var

	def Eg2(self):# mean of g**2
		fun=lambda x:self.std_normal(x)*self.g(self.myTF.TF(x))*self.g(self.myTF.TF(x))
		var,err=integrate.quad(fun,-10,10.)
		return var
