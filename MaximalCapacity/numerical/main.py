import numpy as np
import cPickle as pickle
from classes.transfer_function import * # capaity curves
from classes.learning_rule import * # integrating using einsum
from classes.gaussian_integral import * # integrating using built in scipy
from classes.critical_del0_m import * # capcity curves
import matplotlib.pyplot as plt
#importing data
paramfit = pickle.load(open('parametersFit.p','rb'))

# using the median parameters of the fits
rmax_median = np.median(paramfit[0][0])
beta_median = np.median(paramfit[0][1])
h0_median = np.median(paramfit[0][2])

paramTF = ['sig',1.,beta_median,h0_median] # param TF
tf = TransferFunction(paramTF) 
paramsLR = [26.,26.,0.28,0.28,0,1]
lr = LearningRule(paramsLR,tf)
cr = CriticalOverlap(lr,tf)

# A
if True:
	the_p = np.arange(0,1.1,0.1)
	the_alpha = np.linspace(0.01,.4,200)
	#A - > infty and p free
	cr.MFT_limit(5) 
	the_overlaps = []
	for p in the_p:
		cr.GI.p = p
		del0 = 0.01
		m = 1.
		overlaps= []
		for alpha in the_alpha:
			del0,m = cr.self_consistent(del0,m,alpha)
			overlaps.append(m)
			print 'alpha=',alpha
		print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
		print 'p=',p
		print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
		the_overlaps.append(overlaps)
	pickle.dump((the_p,the_alpha,the_overlaps),open('curves/overlaps_A_inf.p','wr'))


# plot B
if True:
	cr.MFT_limit(3) #A and p free
	the_p = np.logspace(-3,0,400)
	the_A = np.array([5.,6.1,6.95,8.,10.,20,100.,1000.])
	the_capacity = []
	for A in the_A:
		capacity = []
		cr.GI.A_bar = A
		for p in the_p:
			cr.GI.p = p
			capacity.append(cr.capacity())
			print 'A=',cr.GI.A_bar,'p=',cr.GI.p,'alpha_c',capacity[-1]
		the_capacity.append(capacity)

	pickle.dump((the_A,the_p,the_capacity),open('curves/capacity_vs_p_A.p','wr'))


#C
if True:
	# p = 0
	the_A_p_0 = np.logspace(0,3,400)
	the_capacity_p_0 = []

	p = 1e-3
	cr.GI.p = p
	#p = 0  and qf!=1 <--> qf !=  qg
	cr.MFT_limit(4) #A and p free

	#p = 0  and qf!=1 <--> qf !=  qg
	cr.MFT_limit(4) #A, p, qf free
	del_qf = [0.,0.02,0.05,0.1]
	for qf in del_qf:
		cr.GI.qf = 1-p-qf
		capacity = []
		for A in the_A_p_0:
			cr.GI.A_bar = A
			capacity.append(cr.capacity())
			print 'A=',A,'p=',cr.GI.p,'qf=',cr.GI.qf,'alpha_c',capacity[-1]
		the_capacity_p_0.append(capacity)

	pickle.dump((the_A_p_0,del_qf,the_capacity_p_0),open('curves/capacity_vs_A_p0.p','wr'))



#D
if True:
	#value A
	cr.GI.A_bar = 6.95
	the_p = np.logspace(-3,np.log10(0.5),400)
	the_capacity = []
	
	#p = 0  and qf!=1 <--> qf !=  qg
	cr.MFT_limit(4) #A, p and qf free
	del_qf = [0.,0.02,0.05,0.1]#[0.,0.1,0.3,0.5]
	for qf in del_qf:
		capacity = []
		for p in the_p:
			cr.GI.qf = 1-p-qf
			cr.GI.p = p
			capacity.append(cr.capacity())
			print 'A=',cr.GI.A_bar,'qf=',cr.GI.qf,'p=',cr.GI.p,'alpha_c',capacity[-1]
		the_capacity.append(capacity)
	
	pickle.dump((del_qf,the_p,the_capacity),open('curves/capacity_vs_p_qf.p','wr'))

