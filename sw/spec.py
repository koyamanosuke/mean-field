from locale import YESEXPR
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys

texparams = {'ps.useafm': True, 'pdf.use14corefonts': True, 'pdf.fonttype': 42,
             'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath} \usepackage{txfonts} \usepackage{bm}'}
plt.rcParams.update(texparams)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']



from matplotlib.colors import LinearSegmentedColormap
def generate_cmap(colors):
    """自分で定義したカラーマップを返す"""
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v / vmax, c))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def generate_cmapv(colors):
    """自分で定義したカラーマップを返す"""
    values = [v[1] for v in colors]
    vmax = np.max(values)
    vmin = np.min(values)
    color_list = []
    for c in colors:
        color_list.append(((c[1]-vmin)/(vmax-vmin), c[0]))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


cm = generate_cmapv([('black', 0.), ('navy', 0.05), ('#A001B8', 0.1), ('red', 0.2), ('orange', 0.3),
                     ('yellow', 0.6), ('#bbeebb', 0.8), ('#e0ffff', 0.9), ('white', 1)])
win = generate_cmapv([('black', 0.), ('#5C00BE', 0.05), ('#12FFEF', 0.2), ('#9EFF51', 0.3),
                     ('yellow', 0.4), ('#FE4300', 0.5), ('#ED0000', 0.6), ('#B40000', 0.8), ('#7B0000', 1)])



def fl(x, x0, gamma):
    return (1./np.pi)*gamma/((x-x0)*(x-x0) + gamma*gamma)


# gamma = 0.05

ndata = np.loadtxt('data/output/spec.txt')
# damp = np.loadtxt('data/output/damp.txt')

arg = sys.argv
M = int(arg[1])   #  副格子(独立なクラスター)の数
Ne = int(arg[2])  #  励起状態の数
N = M * Ne

# for i in range(1, M * Ne + 1):
#     ndata[:, 11 * M * Ne + i] = ndata[:, 11 * M * Ne + i] - ndata[:, i]

X = ndata[:, 0]
ymax = np.max(ndata[:, M * Ne]) * 1.5
Y = np.linspace(-2, ymax, 4000)
Z = np.empty((Y.size, X.size))

nw = np.empty((N, Y.size, X.size))
nwa = np.empty((N, X.size))
nweight = np.empty((N, Y.size, X.size))

x, y = np.meshgrid(X, Y)

for i in range(N):
    nw[i,:,:], NULL = np.meshgrid(ndata[:,i+1], Y)

for i in range(0, N):
    nwa[i] = (ndata[:, N * 1 + 1 + i] + ndata[:, N * 2 + 1 + i] + ndata[:, N * 3 + 1 + i] + + ndata[:, N * 4 + 1 + i])
    # n = ndata.shape[0]
    # nwa[i] = np.ones(n)
    nweight[i], NULL = np.meshgrid(nwa[i], Y)
    # gamma = damp[:, i] + 0.05
    gamma = 0.1
    gamma, Null = np.meshgrid(gamma, Y)
    Z = Z + nweight[i] * fl(y, nw[i], gamma)

#
fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111)

c = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cm, vmin=0., vmax=10.,
                  linewidth=0, rasterized=True)


r3 = np.sqrt(3)
r2 = np.sqrt(2)
xposi = np.array([-2., 0, 2.])
xposi_l = [r'$-2\pi(1,1,1)$', 0, r'$2\pi(1,1,1)$']
ax.set_xticks(xposi)
ax.set_xticklabels(xposi_l, fontsize=14)
xmax = np.max(ndata[:,0])
ax.set_xlim(xposi[0], xposi[xposi.size-1])


yposi = np.arange(1., ymax, 1.)
ax.set_yticks(yposi)
ax.set_yticklabels(yposi, fontsize=16)
ax.set_ylim(0, ymax)

for x in xposi:
    ax.axvline(x, color='white', ls='-', lw=0.5)
ax.axhline(y=0, color='white', ls='-', lw=0.5) 

# ts_bar = [r"0.0", r"1.0", r"2.0", r"3.0"]
bar = fig.colorbar(c, ax=ax)
# bar.ax.set_yticklabels(ts_bar, fontsize=18)

plt.savefig('fig/spec.pdf', transparent=True, bbox_inches='tight')
