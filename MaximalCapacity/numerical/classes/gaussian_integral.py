from transfer_function import *
from learning_rule import *
import numpy as np
from scipy.stats import multivariate_normal as mvnormal
from scipy import integrate
from scipy.stats import norm
from scipy.special import erf

class GaussianIntegrals:
	'''This is a supper class of gaussian integrals'''

	def __init__(self,LR,TF):
		
		self.myTF=TF # a transfer function class
		self.myLR=LR # learning rule class
		
		# integration variables
		self.dx=0.01 #
		self.x_min=-10.
		self.x_max=10.
		self.xgrid = np.arange(self.x_min,self.x_max,self.dx) #grid on x
		self.normal_pdf =self.myLR.std_normal(self.xgrid) #standard normal pdf
		
		self.deta=0.1 #
		eta_min=0.1
		eta_max=self.myTF.rm	
		self.eta_grid = np.arange(eta_min,eta_max,self.deta) #grid on x
		self.pdf_eta =self.myTF.dist_rates(self.eta_grid) #standard normal pdf
		self.f_eta=self.myLR.f(self.eta_grid)# trasnformtion eta by f
		self.g_eta=self.myLR.g(self.eta_grid)# trasnformtion eta by g
		# defining used vectors/matrices

		# built in integration
		self.nquad=200

		#parameters limit f and g step
		self.p = 0.5
		self.qf = 0.5
		self.A_bar = 6.95

	def update_grids(self,dx,x_min,x_max,deta,Eta_min):
		# integration variables
		self.dx=dx #
		self.x_min=x_min
		self.x_max=x_max
		self.xgrid = np.arange(self.x_min,self.x_max,self.dx) #grid on x
		self.normal_pdf =self.myLR.std_normal(self.xgrid) #standard normal pdf
	
		self.deta=deta #
		eta_min=Eta_min
		eta_max=self.myTF.rm	
		self.eta_grid = np.arange(eta_min,eta_max,self.deta) #grid on x
		self.pdf_eta =self.myTF.dist_rates(self.eta_grid) #standard normal pdf
		self.f_eta=self.myLR.f(self.eta_grid)# trasnformtion eta by f
		self.g_eta=self.myLR.g(self.eta_grid)# trasnformtion eta by g
		# defining used vectors/matrices

	
	# different grid eta and x
	def Int_overlap_v2(self,del0,m):	
		meanfield=self.myLR.Amp*self.f_eta*m
		sdfield=self.myLR.Amp*np.sqrt(del0)*self.xgrid
		steady_state=self.myTF.TF(np.add.outer(meanfield,sdfield))
		g=self.pdf_eta*self.g_eta
		sol=self.dx*self.deta*np.einsum('i,ij,j',g,steady_state,self.normal_pdf)
		return sol

	# integraditon using built in method
	def rbar(self,del0,m,z):
		# First moment for PW linear TF
		sigma=self.myLR.Amp*np.sqrt(del0)
		mu = self.myLR.f(self.myTF.TF(z))*m*self.myLR.Amp
		r=lambda x:self.myLR.std_normal(x)*self.myTF.TF(sigma*x+mu)
		var,err=integrate.fixed_quad(r,-10,10,n=self.nquad)
		return var
	
	def Int_overlap_v3(self,del0,m):	
		overlap = lambda z:self.myLR.std_normal(z)*self.myLR.g(self.myTF.TF(z))*self.rbar(del0,m,z)
		var,err = integrate.quad(overlap,-np.inf,np.inf)
		return var
	
	# limit  step with A, p and qf free
	def Int_overlap_step_LR_qf(self,del0,m):	
		A =  self.A_bar
		p = self.p
		Eg2 = p * (1 - p)
		Ef2 = p * self.qf**2 + (1 - p) * (1 - self.qf)**2
		sigma =  np.sqrt(Eg2/Ef2)
		#print 'simga overlap=',sigma
		mean_field_p = A * self.qf * m * sigma
		mean_field_n = -A * (1-self.qf) * m * sigma
		sd_field = A * np.sqrt(del0) * self.xgrid
		TF_p = self.myTF.TF(mean_field_p + sd_field)
		TF_n = self.myTF.TF(mean_field_n + sd_field)
		sol_p = self.dx * np.einsum('i,i',TF_p,self.normal_pdf)
		sol_n = self.dx * np.einsum('i,i',TF_n,self.normal_pdf)
		return sol_p -  sol_n
	
	# limit A and p free
	def Int_overlap_step_LR(self,del0,m):	
		A =  self.A_bar
		p =  self.p
		mean_field_p = A * (1 - p) * m
		mean_field_n = -A * p * m
		sd_field = A * np.sqrt(del0) * self.xgrid
		TF_p = self.myTF.TF(mean_field_p + sd_field)
		TF_n = self.myTF.TF(mean_field_n + sd_field)
		sol_p = self.dx * np.einsum('i,i',TF_p,self.normal_pdf)
		sol_n = self.dx * np.einsum('i,i',TF_n,self.normal_pdf)
		return sol_p -  sol_n
	
	# derivative m
	def DInt_overlap_step_LR_p0(self,del0,m):	
		A =  self.A_bar
		mean_field = A * m
		sd_field = A * np.sqrt(del0) * self.xgrid
		DTF = self.myTF.DTF(mean_field + sd_field)
		sol = self.dx * np.einsum('i,i',DTF,self.normal_pdf)
		return sol 
	
	# limit A \to infinity and f and g step
	def Int_overlap_large_A(self,del0,m):	
		p =  self.p
		mean_field_p = (1 - p) * m
		mean_field_n = - p * m
		sd_field = np.sqrt(del0) 
		sol_p = 1 - 0.5 * (1 + erf(-mean_field_p/(np.sqrt(2)*sd_field)))
		sol_n = 1 - 0.5 * (1 + erf(-mean_field_n/(np.sqrt(2)*sd_field)))
		return sol_p -  sol_n
	
	# integraditon using built in method
	def r2bar(self,del0,m,z):
		# First moment for PW linear TF
		sigma=self.myLR.Amp*np.sqrt(del0)
		mu = self.myLR.f(self.myTF.TF(z))*m*self.myLR.Amp
		r=lambda x:self.myLR.std_normal(x)*self.myTF.TF(sigma*x+mu)**2
		var,err=integrate.fixed_quad(r,-10,10,n=self.nquad)
		return var
	
	def Int_TF2_v3(self,del0,m):	
		overlap = lambda z:self.myLR.std_normal(z)*self.r2bar(del0,m,z)
		var,err = integrate.quad(overlap,-np.inf,np.inf)
		return var
	
	# limit f and g step A, qf and p free
	def Int_TF2_step_LR_qf(self,del0,m):	
		A =  self.A_bar
		p =  self.p
		Eg2 = p * (1 - p)
		Ef2 = p * self.qf**2 + (1 - p) * (1 - self.qf)**2
		sigma =  np.sqrt(Eg2/Ef2)
		mean_field_p = A * self.qf * sigma
		mean_field_n = -A * (1 - self.qf) * m * sigma
		sd_field = A * np.sqrt(del0) * self.xgrid
		TF2_p = self.myTF.TF(mean_field_p + sd_field)**2
		TF2_n = self.myTF.TF(mean_field_n + sd_field)**2
		sol_p = self.dx*np.einsum('i,i',TF2_p,self.normal_pdf)
		sol_n = self.dx*np.einsum('i,i',TF2_n,self.normal_pdf)
		return p * sol_p + (1 - p) * sol_n 
	
	# Limit f and g step A and p free
	def Int_TF2_step_LR(self,del0,m):	
		A =  self.A_bar
		p =  self.p
		mean_field_p = A * (1 - p) * m
		mean_field_n = -A * p * m
		sd_field = A * np.sqrt(del0) * self.xgrid
		TF2_p = self.myTF.TF(mean_field_p + sd_field)**2
		TF2_n = self.myTF.TF(mean_field_n + sd_field)**2
		sol_p = self.dx*np.einsum('i,i',TF2_p,self.normal_pdf)
		sol_n = self.dx*np.einsum('i,i',TF2_n,self.normal_pdf)
		return p * sol_p + (1 - p) * sol_n 
	
	# limit A \to infinity and f and g step
	def Int_TF2_large_A(self,del0,m):	
		p =  self.p
		mean_field_p = (1 - p) * m
		mean_field_n = - p * m
		sd_field = np.sqrt(del0) 
		sol_p = 1 - 0.5 * (1 + erf(-mean_field_p/(np.sqrt(2)*sd_field)))
		sol_n = 1 - 0.5 * (1 + erf(-mean_field_n/(np.sqrt(2)*sd_field)))
		return p * sol_p + (1 - p) * sol_n 
