import numpy as np
from classes.learning_rule import *
from classes.transfer_function import *
from classes.mft_curves import *
import cPickle as pickle
import matplotlib.pyplot as plt



#importing data
paramfit = pickle.load(open('../parametersFit.p','rb'))

# using the median parameters of the fits
rmax_median = np.median(paramfit[0][0])
beta_median = np.median(paramfit[0][1])
h0_median = np.median(paramfit[0][2])
paramTF = ['sig',rmax_median,beta_median,h0_median] # param TF


tf = TransferFunction(paramTF) 

#learning rule
amp_median = np.median(paramfit[1][0])
qf_median = np.median(paramfit[1][1])
bf_median = np.median(paramfit[1][2])
xf_median = np.median(paramfit[1][3])

# capacity vs bg
if True:
	cap = []
	the_bg = np.logspace(-2,2,200)
	init = 0.01
	for bg in the_bg:
		paramLR = [xf_median,xf_median,bf_median,bg,qf_median,amp_median]
		lr = LearningRule(paramLR,tf)
		curves = MFTCurves(lr,tf)
		alpha_c = curves.capacity_static_MFT(init)
		cap.append(alpha_c)
		print 'bg=',bg

	pickle.dump((the_bg,np.array(cap)),open('curves/bg_vs_cap.p','wb'))


if True:
	cap = []
	the_xg = np.arange(10,40,0.1)
	init = 0.01
	for xg in the_xg:
		paramLR = [xf_median,xg,bf_median,bf_median,qf_median,amp_median]
		lr = LearningRule(paramLR,tf)
		curves = MFTCurves(lr,tf)
		alpha_c = curves.capacity_static_MFT(init)
		cap.append(alpha_c)
		print 'xg=',xg
	pickle.dump((the_xg,np.array(cap)),open('curves/xf_vs_cap.p','wb'))



