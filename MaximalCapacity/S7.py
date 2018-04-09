import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt

#loading
the_A_p,the_p_cap,the_capacity_p = pickle.load(open('numerical/curves/capacity_vs_p_A.p','rb'))
the_A_p_0,del_q_p_0,the_capacity_p_0 = pickle.load(open('numerical/curves/capacity_vs_A_p0.p','rb'))
the_p_overlap,the_alpha,the_overlaps = pickle.load(open('numerical/curves/overlaps_A_inf.p','rb'))
delta_q_qf,the_p_qf,the_capacity_qf = pickle.load(open('numerical/curves/capacity_vs_p_qf.p','rb'))
alpha_max = max(the_capacity_p_0[0])



#plotting
params = {'legend.fontsize': 10,'axes.labelsize': 16,'axes.titlesize':22,\
		'xtick.labelsize':14,'ytick.labelsize':14}
plt.rcParams.update(params)
figure=plt.figure(figsize=(13,10))
figure.subplots_adjust(hspace = 0.45,wspace = 0.3)

fig1=figure.add_subplot(221)

for i in range(0,len(the_overlaps),2):
	fig1.plot(the_alpha,the_overlaps[i],label=r'$p=$'+str(the_p_overlap[i]),lw=2)
fig1.axvline(x=1./np.pi, ymin=0, ymax=1., linewidth=1,ls ='--', color = 'r')
fig1.set_xlim([0,0.42])
fig1.set_xticks([0.2,0.3183,0.4])
fig1.set_xticklabels([0.2,r'$\frac{1}{\pi}$',0.4])
fig1.set_yticks([0.,0.5,1.])
fig1.set_ylim(0,1.2)
fig1.set_xlabel(r'$\alpha$')
fig1.set_ylabel(r'$m_0$')
fig1.legend(loc='upper right')
fig1.set_title('(A)',y=1.04)




# fig b
labels = [r'=5$','=6.1$',r'\approx 7$',r'=8.$',r'=10$',r'=20$',r'=10^{2}$',r'=10^{3}$']
fig2=figure.add_subplot(222)
i=0
for cap in the_capacity_p:
	fig2.semilogx(the_p_cap,cap,label=r'$\bar{A}'+labels[i],lw=2)
	i=i+1
fig2.axhline(y=1./np.pi, xmin=0, xmax=1., linewidth=1,ls ='--', color = 'r')
#fig2.axhline(y=alpha_max, xmin=0, xmax=1., linewidth=1,ls ='--', color = 'r')
fig2.set_ylim([0,1.])
fig2.set_yticks([0.3183,0.4,0.8])
fig2.set_yticklabels([r'$\frac{1}{\pi}$',0.4,0.8])
fig2.set_xlim(0,2.)
fig2.set_xlabel(r'$p$')
fig2.set_ylabel(r'$\alpha_c$')
fig2.legend(loc='upper right')
fig2.set_title('(B)',y=1.04)

# fig c
sigma = lambda p,qf:(p * (1-p))/(p*qf**2+(1-p)*(1-qf)**2)

fig3=figure.add_subplot(223)
p = 1e-3
ylabels= [r'$\frac{1}{\pi}$']
yticks= [ ]
for l in range(len(the_capacity_p_0)-1):
	fig3.semilogx(the_A_p_0,the_capacity_p_0[l],label =r'$\Delta_q =$'+str(del_q_p_0[l]) )
	qf = 1 - p - del_q_p_0[l]
	sigma2 = sigma(p,qf)
	fig3.axhline(y=sigma(p,qf)/np.pi, xmin=0, xmax=2000, linewidth=1,ls ='--', color = 'r')
	yticks.append(sigma2/np.pi)
	if 0<l:
		ylabels.append(r'$\frac{'+str(round(sigma2,1))+'}{\pi}$')
ylabels.append(0.4)
yticks.append(0.4)
ylabels.append(0.8)
yticks.append(0.8)
fig3.axhline(y=1./np.pi, xmin=0, xmax=2000, linewidth=1,ls ='--', color = 'r')
#fig3.axhline(y=alpha_max, xmin=0, xmax=1., linewidth=1,ls ='--', color = 'r')
fig3.set_xlabel(r'$\bar{A}$')
fig3.set_ylim([0,1.])
fig3.set_xlim([0,1000])
fig3.set_yticks(yticks)
fig3.set_yticklabels(ylabels)
fig3.set_ylabel(r'$\alpha_c$')
fig3.legend(loc='upper right')
fig3.set_title('(C)',y=1.04)


# fig d
fig4=figure.add_subplot(224)
i = 0
for capacity_qf in the_capacity_qf[0:3]:
	fig4.semilogx(the_p_qf,capacity_qf,label = r'$\Delta_q = $'+str(delta_q_qf[i]))
	i+=1
#fig4.axhline(y=1./np.pi, xmin=0, xmax=1., linewidth=1,ls ='--', color = 'r')
#fig4.axhline(y=alpha_max, xmin=0, xmax=1., linewidth=1,ls ='--', color = 'r')
fig4.set_ylim([0,1.])
fig4.set_yticks([0.4,0.8])
fig4.set_yticklabels([0.4,0.8])
fig4.set_xlim(0,1.)
fig4.set_xlabel('$p$')
fig4.set_ylabel(r'$\alpha_c$')
fig4.legend(loc='upper right')
fig4.set_title('(D)',y=1.04)

plt.savefig('S7.pdf', bbox_inches='tight')
