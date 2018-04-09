import numpy as np
import cPickle as pickle
from cycler import cycler
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

''' This code calculate the maxcapacity for smooth functions 
of f and g'''

#
overlaps_N = []
PATH='numeric/mft_sim/overlaps/overlaps_sim'
all_m_sim = []
all_M_sim = []

n_real = 4

num_patterns,ord_params = pickle.load(open(PATH+'.p','rb'))
for i in range(0,n_real):
	all_m_sim.append(ord_params[i][0])
	all_M_sim.append(ord_params[i][1])


mean_m = np.mean(all_m_sim,axis=0)
sd_m = np.std(all_m_sim,axis=0)

mean_M = np.mean(all_M_sim,axis=0)
sd_M = np.std(all_M_sim,axis=0)

myalp=(1./250.) * num_patterns

# loading theoretical overlap curve
PATH='numeric/mft_curves/curves/'
thealpha,m_theory,M_theory=pickle.load(open(PATH+'overlap_mft.p','rb'))

#loading dynamics overlaps and  neurons
PATH2='numeric/network_simulation/overlaps/'
mydyn=[]
mym=[]
replicas=10
mydyn = pickle.load(open(PATH2+'the_dynamics.p','rb'))
mym = pickle.load(open(PATH2+'the_overlaps.p','rb'))
mysoldyn_fam=mydyn[0]
time_fam=np.linspace(0,8.5,len(mysoldyn_fam[1000:-1,0]))

fig = plt.figure(figsize=(24, 18))
gs = gridspec.GridSpec(10,10)
gs.update(wspace=0.05,hspace=0.05)
for i in range(10):
	for j in range(10):
		ax=plt.subplot(gs[i,j])
		l=10*i+j
		print l
		ax.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
		ax.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
		for k in range(replicas):
			ax.plot(time_fam,mydyn[k][1000:-1,l],color='r',lw=2,alpha=0.1)
		average_trace=np.mean(np.array([mydyn[r][1000:-1,l] for r in range(len(mydyn))]),axis=0)
		ax.plot(time_fam,average_trace,color='k',lw=2,alpha=1.)
		ax.set_xlim([0,4])
		ax.set_ylim([0,80])
		ax.set_yticks([])
		ax.set_xticks([])
		#fg3.set_ylabel(r'Rate (Hz)',fontsize=25)
		#fg3.set_xlabel(r'Time (s)',fontsize=25)
		#fg3.tick_params(labelsize=20)

plt.savefig('fig_all_neurons.pdf', bbox_inches='tight')
plt.close()

#--------------------------------------------------------------------
#--------------------Main Figure ------------------------------------
#--------------------------------------------------------------------



rep = replicas # number of replicas simulation
initdist=plt.figure(figsize=(14,7))
initdist.subplots_adjust(wspace=.4) # vertical space bw figures
initdist.subplots_adjust(hspace=.7) # vertical space bw figures



myind=0#index of realization
num_neu=20
colormap = plt.cm.jet
plt.rc('axes', prop_cycle=(cycler('color', [colormap(i) for i in np.linspace(0.1,0.9,20)])))
fg3=initdist.add_subplot(231)
fg3.plot(time_fam,mydyn[myind][1000:-1,0:num_neu],lw=2,alpha=0.5)
fg3.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.set_xlim([0,5.5])
fg3.set_ylim([0,80])
fg3.set_yticks([40,80])
fg3.set_xticks([0,2.5,5])
#fg1.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
fg3.set_ylabel(r'Rate (Hz)',fontsize=25)
fg3.set_xlabel(r'Time (s)',fontsize=25)
fg3.tick_params(labelsize=20)
fg3.set_title('(A)',fontsize=25,y=1.06)

fg3=initdist.add_subplot(232)
fg3.plot(time_fam,mym[myind][1000:-1,0],color='b',lw=5,alpha=1.,label=r'Sim. $m_1$')
fg3.plot(time_fam,mym[myind][1000:-1,1],color='g',lw=5,alpha=0.4,label=r'Sim. $m_k$')
fg3.plot(time_fam,mym[myind][1000:-1,2:-1],color='g',lw=5,alpha=0.4)
#fg3.axhline(y=overlapFamDelayTheo,xmin=5./11.,xmax=1.,color='r',lw=6,alpha=0.5,label=r'MFT')
#fg3.axhline(y=overlapFamPresTheo,xmin=4./11,xmax=5./11.,color='r',lw=6,alpha=0.5)
fg3.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.set_ylim([-0.1,1.2])
fg3.set_yticks([0,0.5,1])
fg3.set_xlim([0,5.5])
fg3.set_xticks([0,2.5,5])
#fg1.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
fg3.set_ylabel(r'Overlap '+r'($m$)',fontsize=25)
fg3.set_xlabel(r'Time (s)',fontsize=25)
fg3.tick_params(labelsize=20)
fg3.legend(loc='upper right',numpoints=1,prop={'size':13.5})
fg3.set_title('(B)',fontsize=25,y=1.06)

#overlap
#m_theory=np.array(m_theory)
fg2=initdist.add_subplot(233)
fg2.plot(thealpha[m_theory>0.01],m_theory[m_theory>0.01],lw=10,ls='-',color='k',alpha=0.5,label='MFT')
fg2.axvline(x=thealpha[m_theory>0.01][-1],ymin=0,ymax=m_theory[m_theory>0.01][-1]/1.2,lw=3,color='k',ls=':')
fg2.errorbar(myalp[0:8],mean_m[0:8],yerr=sd_m[0:8],alpha=0.5,fmt='o',ms=15,markerfacecolor='r',markeredgecolor='r',ecolor='r',label=r'Fixed point')
fg2.errorbar(myalp[8:20],mean_m[8:20],yerr=sd_m[8:20],alpha=0.5,fmt='s',ms=15,markerfacecolor='b',markeredgecolor='b',ecolor='b',label=r'Chaos')
fg2.set_xlim([0,0.9])
fg2.set_ylim([0,1.2])
fg2.set_yticks([0.5,1])
fg2.set_xticks([0,0.2,0.4,0.6,0.8])
fg2.set_xlabel(r'Memory load ($\alpha$)',fontsize=25)
fg2.set_ylabel(r'Overlap ($m$)',fontsize=25)
fg2.tick_params(labelsize=20)
fg2.set_title('(C)',fontsize=25,y=1.06)
leg = fg2.legend(loc='upper right',numpoints=1,prop={'size':10.5},markerscale=0.8)
leg_lines = leg.get_lines()
plt.setp(leg_lines, linewidth=6)
#lgnd.legendHandles[2]._legmarker.set_markersize(8)


colormap = plt.cm.jet
fg3=initdist.add_subplot(234)
index_neuron_d=3
fg3.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
#plt.rc('axes', prop_cycle=(cycler('color', [colormap(i) for i in np.linspace(0.1,0.9,3)])))
for i in range(replicas):
	fg3.plot(time_fam,mydyn[i][1000:-1,index_neuron_d],color='b',lw=2,alpha=0.1)
plt.rc('axes', prop_cycle=(cycler('color', [colormap(i) for i in np.linspace(0.1,0.9,3)])))
average_trace_d=np.mean(np.array([mydyn[r][1000:-1,index_neuron_d] for r in range(rep)]),axis=0)

fg3.plot(time_fam,average_trace_d,color='b',lw=2,alpha=1.)
fg3.set_xlim([0,4])
fg3.set_ylim([0,80])
fg3.set_yticks([40,80])
fg3.set_xticks([0,2,4])
fg3.set_ylabel(r'Rate (Hz)',fontsize=25)
fg3.set_xlabel(r'Time (s)',fontsize=25)
fg3.tick_params(labelsize=20)
fg3.set_title('(D)',fontsize=25,y=1.06)

colormap = plt.cm.jet
fg3=initdist.add_subplot(235)
index_neuron_e=40
fg3.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
for i in range(replicas):
	fg3.plot(time_fam,mydyn[i][1000:-1,index_neuron_e],color='r',lw=2,alpha=0.1)
average_trace_e=np.mean(np.array([mydyn[r][1000:-1,index_neuron_e] for r in range(rep)]),axis=0)
fg3.plot(time_fam,average_trace_e,color='r',lw=2,alpha=1.)
fg3.set_xlim([0,4])
fg3.set_ylim([0,80])
fg3.set_yticks([40,80])
fg3.set_xticks([0,2.,4.])
fg3.set_ylabel(r'Rate (Hz)',fontsize=25)
fg3.set_xlabel(r'Time (s)',fontsize=25)
fg3.tick_params(labelsize=20)
fg3.set_title('(E)',fontsize=25,y=1.06)

colormap = plt.cm.jet
fg3=initdist.add_subplot(236)
index_neuron_f=55#23#89
fg3.axvline(x=2.,ymin=0,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
fg3.axvline(x=2.5,ymin=0.,ymax=5.,color='gray',lw=3,alpha=1.,ls='--')
for i in range(replicas):
	fg3.plot(time_fam,mydyn[i][1000:-1,index_neuron_f],color='g',lw=2,alpha=0.1)
average_trace_f=np.mean(np.array([mydyn[r][1000:-1,index_neuron_f] for r in range(rep)]),axis=0)
fg3.plot(time_fam,average_trace_f,color='g',lw=2,alpha=1.)
fg3.set_xlim([0,4.])
fg3.set_ylim([0,80])
fg3.set_yticks([40,80])
fg3.set_xticks([0,2,4])
fg3.set_ylabel(r'Rate (Hz)',fontsize=25)
fg3.set_xlabel(r'Time (s)',fontsize=25)
fg3.tick_params(labelsize=20)
fg3.set_title('(F)',fontsize=25,y=1.06)


plt.savefig('fig.pdf', bbox_inches='tight')
plt.close()


