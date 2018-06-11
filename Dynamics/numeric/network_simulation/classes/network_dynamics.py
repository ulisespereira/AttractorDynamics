import numpy as np
from scipy import sparse 

class NetworkDynamics:
	'''This class creates the connectivity matrix'''

	def __init__(self,LR,TF,connectivity,patterns_current):
		
		#tranfer function and learning rule
		self.myLR=LR
		self.myTF=TF

		#
		self.dt=0.5 # dt integration
		self.tau=20. # 20ms 

		#input current
		self.Input=0.
		self.connectivity = connectivity
		self.patterns_fr = self.myTF.TF(patterns_current)
		self.patterns_current = patterns_current	
		self.amp_stim = 1.#5
	def fieldDynamics(self,u,t):
		return (1./self.tau)*(-u+self.myTF.TF(self.Input+self.connectivity.dot(u)))#-1.*np.mean(u)))
	# dynamics 
	def dynamics(self,period,u_init,patterns_fr):
		T=period
		n_neurons=100#n neurons to save
		un=u_init#self.myTF.TF(u_init) 
		p,N=patterns_fr.shape
		
		mysol=[] #neurons dynammics
		q_ord_param=[] # overlap
		m_ord_param=[] # normalized overlap
		
		
		mysol.append(un[0:n_neurons])
		overlap=np.array([np.mean(np.multiply(self.myLR.g(patterns_fr[i,:]),un)) for i in range(p)])
	
		q_ord_param.append(overlap)
		m_ord_param.append(overlap/(np.sqrt(self.myLR.intg2)*np.std(un)))
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
			q_ord_param.append(overlap)
			m_ord_param.append(overlap/(np.sqrt(self.myLR.intg2)*np.std(un)))
			print 't',t,'of',T	
		
		return np.array(q_ord_param),np.array(m_ord_param),np.array(mysol),un
	
	# delay response task
	def DMS(self,T_bckg,T_pres,T_delay,u_init):
		#background period
		print 'Background period'
		self.Input = 0.
		q_bckg,m_bckg,sol_bckg,un = self.dynamics(T_bckg,u_init,self.patterns_fr)
		#presentation period
		print 'Presentation period'
		self.Input = self.amp_stim * self.patterns_current[0]
		q_pres,m_pres,sol_pres,un = self.dynamics(T_pres,un,self.patterns_fr)
		#delay period
		print 'Delay period'
		self.Input = 0.
		q_delay,m_delay,sol_delay,un = self.dynamics(T_delay,un,self.patterns_fr)
		#concatenate
		q = np.concatenate((q_bckg,q_pres,q_delay),axis = 0)
		m = np.concatenate((m_bckg,m_pres,m_delay),axis = 0)
		sol = np.concatenate((sol_bckg,sol_pres,sol_delay),axis = 0)
		
		return q,m,sol
