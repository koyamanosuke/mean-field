import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.colors as matcol
from matplotlib.cm import ScalarMappable, get_cmap

texparams = {'ps.useafm': True, 'pdf.use14corefonts': True, 'pdf.fonttype': 42,
             'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath} \usepackage{txfonts} \usepackage{bm}'}
plt.rcParams.update(texparams)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

arg = sys.argv
M = int(arg[1])
N = int(arg[2])  #局所励起状態の数

w1 = np.loadtxt('spec.txt')

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(1, M * N):
    ax.plot(w1[:, 0],w1[:, i], color='orange', lw=1.)

ax.plot(w1[:,0],w1[:, M * N],label=r'$SWenergy$', color='orange', lw=1.)

ax.legend(fontsize=20)




#plt.xlabel('kx')
ax.set_ylabel(r'$\omega$ ', fontsize=27)
ax.set_xticks([0., 4./3., 5./3., 8./3.])
ax.set_xticklabels([r"$\Gamma$",r"$K$",r"$M$",r"$\Gamma$"],fontsize=27)

xmax = np.max(w1[:,0])
ax.set_xlim(0, xmax)
ymax = np.max(w1[:, M * N]) + 0.5
ax.set_ylim(0, ymax)

yposi = np.arange(0, ymax, 1)
ax.set_yticks(yposi)
ax.set_yticklabels(yposi, fontsize=20)


ax.vlines(0, 0, ymax, linestyle='dashed', linewidth=0.3,) 
ax.vlines(4./3., 0, ymax, linestyle='dashed', linewidth=0.3) 
ax.vlines(8./3., 0, ymax, linestyle='dashed', linewidth=0.3)
# ax.vlines(1, 0, 4, linestyle='dashed',color='red', linewidth=0.8)
# ax.vlines(2.1666, 0, 4, linestyle='dashed', color='green', linewidth=0.8)
 

# plt.yticks(color="None")

ax.tick_params(axis='both', which='both', direction='in',bottom=True, top=True, left=True, right=True)



plt.savefig('khsw.pdf',transparent=True,bbox_inches='tight')