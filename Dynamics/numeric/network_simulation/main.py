import numpy as np
from classes.learning_rule import *
from classes.connectivity import *
from classes.transfer_function import *
from classes.network_dynamics import *
import cPickle as pickle
import multiprocessing as mt
import matplotlib.pyplot as plt
import os
directory = 'overlaps'
if not os.path.exists(directory):
	    os.makedirs(directory)

# fixed-point or chaotic attractors as retrival state
TypeDynamics = 'chaos'
#TypeDynamics = 'fixedpoint'

#importing parameters data
paramfit = pickle.load(open('../parametersFit.p','rb'))

# using the median parameters of the fits
rmax_median = np.median(paramfit[0][0])
beta_median = np.median(paramfit[0][1])
h0_median = np.median(paramfit[0][2])
paramTF = ['sig',rmax_median,beta_median,h0_median] # param TF


# transfer function
tf = TransferFunction(paramTF) 



amp_median = np.median(paramfit[1][0])
qf = np.median(paramfit[1][1])#0.65
bf_median = np.median(paramfit[1][2])
xf = np.median(paramfit[1][3])#22.
paramLR = [xf,xf,bf_median,bf_median,qf,amp_median]  #learning rule

if  TypeDynamics =='chaos':
	amp_median = 3 * amp_median
	p =int(0.56 * 250)#70-80
else:
	amp_median = amp_median	
	p = 30


lr = LearningRule(paramLR,tf)

print 'Parameter Values:'
print 'qg=',lr.qg
print 'xg=',lr.xg
print 'bg=',lr.betag
print 'A=',lr.Amp
print 'qf=',lr.qf
print 'xf=',lr.xf
print 'bf=',lr.betaf



# number of realizations
#connectivity
paramSim = [50000,0.005,p] #N,c,p
random_seed = 7 # random seed connectivity
conn = ConnectivityMatrix(lr,tf,paramSim,random_seed)
matrix = conn.connectivity_generalized_hebbian()
# dynamics
patterns_current = conn.patterns_current


the_overlaps = []
the_dynamics = []
n_real = 10 # number of realizations
i=0
while i<=n_real:
#for i in range(100):
	dyn = NetworkDynamics(lr,tf,matrix,patterns_current)
	u_init = np.random.normal(0,1,paramSim[0])
	q,m,sol = dyn.DMS(2500,500,6000,u_init)
	#u_init = patterns_current[0]
	#q,m,sol = dyn.DMS_Short(500,3000,u_init)
	print 'Value end m:',m[-1,0]
	print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	print 'The value of i=',i
	print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	if 0.2<m[-1,0]:
		the_dynamics.append(sol)
		pickle.dump(the_dynamics,open('overlaps/the_dynamics.p','wb'))
		i+=1


the_overlaps.append(m)
pickle.dump(the_overlaps,open('overlaps/the_overlaps.p','wb'))


