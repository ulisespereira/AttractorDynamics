import numpy as np
from scipy import sparse 

class NetworkDynamics:
	'''This class creates the connectivity matrix'''

	def __init__(self,LR,TF,connectivity):
		
		#tranfer function and learning rule
		self.myLR=LR
		self.myTF=TF

		#
		self.dt=0.5 # dt integration
		self.tau=20. # 20ms 

		#input current
		self.Input=0.
		self.connectivity = connectivity
	
	def fieldDynamics(self,u,t):
		return (1./self.tau)*(-u+self.myTF.TF(self.Input+self.connectivity.dot(u)))#-1.*np.mean(u)))
	# dynamics 
	def dynamics(self,period,u_init,patterns_fr):
		T=period
		n_neurons=1#n neurons to save
		un=self.myTF.TF(4*u_init) #initial condition
		p,N=patterns_fr.shape
		
		mysol=[] #neurons dynammics
		q_ord_param=[] # overlap
		m_ord_param=[] # normalized overlap
		M_ord_param=[] # normalized overlap
		
		
		mysol.append(un[0:n_neurons])
		overlap=np.array([np.mean(np.multiply(self.myLR.g(patterns_fr[i,:]),un)) for i in range(p)])
		r2 = np.sum(un**2)/float(N)
	
		q_ord_param.append(overlap)
		m_ord_param.append(overlap/(np.sqrt(self.myLR.intg2)*np.std(un)))
		M_ord_param.append(r2)
		
		t=0
		while t<=T:
			#if t<T/3:
			##	self.Input = u_init
			#else:
			#	self.Input = 0.
			un=un+self.dt*self.fieldDynamics(un,t)
			t=t+self.dt
			mysol.append(un[0:n_neurons])
			
			#overlap with all the other patterns
			overlap=np.array([np.mean(np.multiply(self.myLR.g(patterns_fr[i,:]),un)) for i in range(p)])
			r2 = np.sum(un**2)/float(N)
			q_ord_param.append(overlap)
			m_ord_param.append(overlap/(np.sqrt(self.myLR.intg2)*np.std(un)))
			M_ord_param.append(r2)
			
			print t,' of ',T
		
		q_ord_param = np.array(q_ord_param)
		m_ord_param = np.array(m_ord_param)
		mysol = np.array(mysol)
		M_ord_param = np.array(M_ord_param)
		return q_ord_param,m_ord_param,M_ord_param,mysol
