import numpy as np
import matplotlib.pyplot as plt
import sys
# import matplotlib.colors as matcol
# from matplotlib.cm import ScalarMappable, get_cmap

texparams = {'ps.useafm': True, 'pdf.use14corefonts': True, 'pdf.fonttype': 42,
             'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath} \usepackage{txfonts} \usepackage{bm}'}
plt.rcParams.update(texparams)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

arg = sys.argv
M = int(arg[1])   #副格子
Ne = int(arg[2])  #局所励起状態の数

w1 = np.loadtxt('data/output/spec.txt')

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)

for i in range(1, M * Ne):
    ax.plot(w1[:, 0], w1[:, i], color='orange', ls="solid", lw=1.)

ax.plot(w1[:, 0], w1[:, M * Ne],label=r'$\omega_{LSW}$', color='orange', ls="solid", lw=1.)



ax.legend(fontsize=16, ncol=3, loc="upper left", frameon=False)

#plt.xlabel('kx')
ax.set_ylabel(r'$\omega$ ', fontsize=18)

r3 = np.sqrt(3)
r2 = np.sqrt(2)
xposi = np.array([-2., 0, 2.])
xposi_l = [r'$-2\pi(1,1,1)$', 0, r'$2\pi(1,1,1)$']
ax.set_xticks(xposi)
ax.set_xticklabels(xposi_l, fontsize=14)
xmax = np.max(w1[:,0])
ax.set_xlim(xposi[0], xposi[xposi.size-1])


ymax = np.max(w1[:, M * Ne]) * 1.5
yposi = np.arange(1., ymax, 1.)
ax.set_yticks(yposi)
ax.set_yticklabels(yposi, fontsize=16)
ax.set_ylim(0, ymax)




for x in xposi:
    ax.axvline(x, color='black', ls='--', lw=0.5)
ax.axhline(y=0, color='black', ls='--', lw=0.5) 
# ax.vlines(1, 0, 4, linestyle='dashed',color='red', linewidth=0.8)
# ax.vlines(2.1666, 0, 4, linestyle='dashed', color='green', linewidth=0.8)
 

# plt.xticks(color="None")
# plt.yticks(color="None")

ax.tick_params(axis='both', which='both', direction='in',bottom=True, top=True, left=True, right=True)



plt.savefig('fig/sw.pdf',transparent=True,bbox_inches='tight')