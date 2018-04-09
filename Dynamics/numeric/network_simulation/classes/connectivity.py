import numpy as np
from scipy.stats import bernoulli
from scipy import sparse 
import time

class ConnectivityMatrix:
	'''This class creates the connectivity matrix'''

	def __init__(self,LR,TF,param_net,random_seed):
		np.random.seed(random_seed) # fixed the seed 
		
		#tranfer function and learning rule
		self.myTF=TF
		self.myLR=LR

		# parameters for the dynamics
		self.N=int(param_net[0])
		self.c=param_net[1]
		self.p=int(param_net[2])			
		
		self.patterns_current = np.random.normal(0.,1., size=(self.p,self.N))
		self.patterns_fr = self.myTF.TF(self.patterns_current)
	
	def seed(self,semilla):
		np.random.seed(semilla)
		self.patterns_current = np.random.normal(0.,1., size=(self.p,self.N))
		self.patterns_fr = self.myTF.TF(self.patterns_current)

	def connectivity_generalized_hebbian(self):

		patterns_pre=self.myLR.g(self.patterns_fr)
		patterns_post=self.myLR.f(self.patterns_fr)	

		print 'Patterns created. N patterns:',self.p
		#number of entries different than zero
		#N2bar=np.random.binomial(self.N*self.N,self.c)
		#row_ind=np.random.randint(0,high=self.N,size=N2bar)
		#column_ind=np.random.randint(0,high=self.N,size=N2bar)
		
		rv=bernoulli(1).rvs
		indexes=sparse.find(sparse.random(self.N,self.N,density=self.c,data_rvs=rv))
		
		row_ind=indexes[0]
		column_ind=indexes[1]
		N2bar = len(indexes[1])
		print 'Structural connectivity created'
		
		dN=300000
		n=N2bar/dN
		connectivity=np.array([])
		for l in range(n):
			# fast way to write down the outer product learning
			con_chunk=np.einsum('ij,ij->j',patterns_post[:,row_ind[l*dN:(l+1)*dN]],patterns_pre[:,column_ind[l*dN:(l+1)*dN]])
			connectivity=np.concatenate((connectivity,con_chunk),axis=0)
			print 'Synaptic weights created:',100.*(l)/float(n),'%'
		con_chunk=np.einsum('ij,ij->j',patterns_post[:,row_ind[n*dN:N2bar]],patterns_pre[:,column_ind[n*dN:N2bar]])
		print 'Synaptic weights created:',100.,'%'
		connectivity=np.concatenate((connectivity,con_chunk),axis=0)		
		connectivity=(self.myLR.Amp/(self.c*self.N))*connectivity
		print 'Synaptic weights created'

		connectivity=sparse.csr_matrix((connectivity,(row_ind,column_ind)),shape=(self.N,self.N))
		print 'connectivity created'

		return connectivity
	


		
	


