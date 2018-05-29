import numpy as np
from mpl_toolkits.mplot3d import axes3d
import cPickle as pickle
from matplotlib import ticker
import matplotlib as mpl
import matplotlib.pyplot as plt


#importing data
paramfit = pickle.load(open('numerical/parametersFit.p','rb'))


#learning rule
bf_median = np.median(paramfit[1][2])
xf_median = np.median(paramfit[1][3])



# plot settings
plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['lines.dashed_pattern'] = [6, 6]
mpl.rcParams['lines.dashdot_pattern'] = [3, 5, 1, 5]
mpl.rcParams['lines.dotted_pattern'] = [1, 3]
mpl.rcParams['lines.scale_dashes'] = False


# mft curves
PATH = 'numerical/mft_curves/curves/'
the_bg,cap_bg = pickle.load(open(PATH+'bg_vs_cap.p','rb'))
the_xf,cap_xf = pickle.load(open(PATH+'xf_vs_cap.p','rb'))


#---------------------------------------------------------------------------------------------------
#------------------Plotting ------------------------------------------------------------------------


initdist=plt.figure(figsize=(21,7))
initdist.subplots_adjust(wspace=.4) # vertical space bw figures
initdist.subplots_adjust(hspace=.3) # vertical space bw figures



#plot B
fg2=initdist.add_subplot(121)
fg2.semilogx(the_bg,cap_bg,lw=10,color='b',alpha=0.7)
fg2.axvline(x=bf_median,ymin=0,ymax=1,lw=7,color='r',ls='--')
fg2.set_xlim([0.01,100])
fg2.set_ylim([0,0.65])
fg2.set_yticks([0.2,0.4,0.6])
fg2.set_xticks([0.01,0.1,bf_median,1,10,100])
fg2.set_xticklabels([r'$10^{-2}$',r'$10^{-1}$',r'$\beta_f$',r'$10^{0}$',r'$10^1$',r'$10^2$'])
fg2.set_xlabel(r'Slope $\beta_g$ (s)',fontsize=40)
fg2.set_ylabel(r'Capacity ($\alpha_c$)',fontsize=40)
fg2.tick_params(labelsize=35)
fg2.set_title('(A)',fontsize=50,y=1.06)
fg2.legend(loc=1,numpoints=1,prop={'size':22})

#plot C
fg3=initdist.add_subplot(122)
fg3.plot(the_xf,cap_xf,lw=10,color='g',alpha=0.7)
fg3.axvline(x=xf_median,ymin=0,ymax=1,lw=7,color='r',ls='--')
fg3.set_xlim([10,32])
fg3.set_ylim([0,0.65])
fg3.set_yticks([0.2,0.4,0.6])
fg3.set_xticks([10,15,20,25,xf_median,30])
fg3.set_xticklabels([r'$10$',r'$15$',r'$20$',r'$25$',r'$x_f$',r'$30$'])
fg3.set_xlabel(r'Threshold $x_g$ (Hz)',fontsize=40)
fg3.set_ylabel(r'Capacity ($\alpha_c$)',fontsize=40)
fg3.tick_params(labelsize=35)
fg3.set_title('(B)',fontsize=50,y=1.06)
fg3.legend(loc=1,numpoints=1,prop={'size':22})

plt.savefig('fig4.pdf', bbox_inches='tight')

