import numpy as np
from classes.learning_rule import *
from classes.transfer_function import *
from classes.mft_curves import *
import cPickle as pickle
import matplotlib.pyplot as plt

''' overlap curve '''


#importing data
paramfit = pickle.load(open('../parametersFit.p','rb'))

# using the median parameters of the fits
rmax_median = np.median(paramfit[0][0])
beta_median = np.median(paramfit[0][1])
h0_median = np.median(paramfit[0][2])
paramTF = ['sig',rmax_median,beta_median,h0_median] # param TF


tf = TransferFunction(paramTF) 

#learning rule
amp_median = 3 * np.median(paramfit[1][0])
qf = np.median(paramfit[1][1])#0.65
bf_median = np.median(paramfit[1][2])
xf = np.median(paramfit[1][3])#22.

# overlap curve 
if True:
	paramLR = [xf,xf,bf_median,bf_median,qf,amp_median]  #learning rule
	lr = LearningRule(paramLR,tf)
	curves = MFTCurves(lr,tf)
	#alpha_c = curves.capacity_static_MFT(0.01)
	the_alpha = np.arange(0.005,0.7,0.005)
	overlap = []
	r2 = []
	q = 10. #initial q
	for alpha in the_alpha:
		the_ov,del0 = curves.overlap_static_MFT(alpha,q_init=q)
		m,q = the_ov
		overlap.append(m)
		r2.append(del0/(alpha * curves.GI.myLR.gamma))
	pickle.dump((the_alpha,np.array(overlap),np.array(r2)),open('curves/overlap_mft.p','wb'))
