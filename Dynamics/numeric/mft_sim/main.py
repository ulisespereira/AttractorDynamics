import numpy as np
from classes.learning_rule import *
from classes.connectivity import *
from classes.transfer_function import *
from classes.network_dynamics import *
import cPickle as pickle
import multiprocessing as mt
import matplotlib.pyplot as plt
''' Generating 100 realizations of the overlap curve '''


#importing data
paramfit = pickle.load(open('../parametersFit.p','rb'))

# using the median parameters of the fits
rmax_median = np.median(paramfit[0][0])
beta_median = np.median(paramfit[0][1])
h0_median = np.median(paramfit[0][2])
paramTF = ['sig',rmax_median,beta_median,h0_median] # param TF

tf = TransferFunction(paramTF) 

#learning rule
amp_median = 3. * np.median(paramfit[1][0])
qf_median = np.median(paramfit[1][1])
bf_median = np.median(paramfit[1][2])
xf_median = np.median(paramfit[1][3])
paramLR = [xf_median,xf_median,bf_median,bf_median,qf_median,amp_median]  #learning rule

lr = LearningRule(paramLR,tf)


#period
period = 1400
# simulations
n_cores = 4 # put here number cores your computer
# simulating
N = 50000
c = 0.005
K = N * c
num_patterns = int((K/250.))*np.arange(10,210,10)
n_realizations = 51
ord_params = []
for i in range(1,n_realizations+1,n_cores):
	print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	print 'Realization:',i,' of ',n_realizations
	print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	def sim(x):
		overlaps = []
		r2 = []
		np.random.seed()
		for p in num_patterns:
			#connectivity
			paramSim = [N,c,p] #N,c,p
			conn = ConnectivityMatrix(lr,tf,paramSim)
			matrix = conn.connectivity_generalized_hebbian()
			# dynamics
			dyn = NetworkDynamics(lr,tf,matrix)
			u_init = conn.patterns_current[0]
			patterns_fr = np.array([conn.patterns_fr[0]])
			q,m,M,sol = dyn.dynamics(period,u_init,patterns_fr)
			overlaps.append(np.mean(m[500:-1][0]))
			r2.append(np.mean(M[500:-1]))
			print 'overlap values:'
			print overlaps
			print 'M values:'
			print r2
		return (overlaps,r2)
	Input = [s for s in range(n_cores)]
	pool = mt.Pool(processes=n_cores)
	the_overlaps =  pool.map(sim,Input)
	for l in range(n_cores):
		ord_params.append((the_overlaps[l][0],the_overlaps[l][1]))		
pickle.dump((num_patterns,ord_params),open('overlaps/overlaps_sim.p','wb'))
