import matplotlib.pyplot as plt
import numpy as np
import brewer2mpl
import h5py as h5

#set2 = brewer2mpl.get_map('Paired', 'qualitative', 12).mpl_colors
set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
set2.insert(0,'crimson')
set2.append('dodgerblue')

cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360,3460,3560]
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']

h = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/1/traub_syn.h5', 'r')
inh_cells = ['bask23','axax23','LTS23', 'bask56', 'axax56', 'LTS56',
             'nRT']

cells = np.array(cell_range)/ 10
cells_spike = np.diff(cells)
fig, ax = plt.subplots(1)

intra_cell_gap = 6
intra_pop_gap = 15
count = 0
for ii,pop_name in enumerate(pop_names):
    cell_count = cells_spike[ii]
    all_spikes = h['/data/event/'+pop_name+'/spikes'].value
    int_spikes = all_spikes[0:cell_count]
    color = set2[ii]
    if pop_name in inh_cells:
        marker = 'v'
    else:
        marker = '^'
    for kk in range(cell_count):
        x = int_spikes[kk]
        y = np.zeros_like(x) + count
        ax.scatter(x, y, s=7, alpha=0.8,facecolor=color,
                   marker=marker, linewidth='0.2', edgecolor='gray')
        count -= intra_cell_gap
    count -= intra_pop_gap
mid_pts = cells_spike / 2
tick_pts = -1*(mid_pts + cells[:-1])*intra_cell_gap
for ii in range(len(pop_names)):
    tick_pts[ii] -= (intra_pop_gap*ii) 
ax.set_yticks(tick_pts)
ax.set_yticklabels(pop_names)

h_curr = h['/data/uniform/'+pop_name+'/i']
tend = h_curr.shape[1] * float(h_curr.attrs.get('dt'))
tstart = float(h_curr.attrs.get('tstart'))
tunit = h_curr.attrs.get('tunit')
ax.set_xticks(np.linspace(tstart,tend,10))
ax.set_xlabel(tunit)
plt.ylim((-2400 ,50))
plt.title('Spikes', fontweight="bold", fontsize=12)
plt.show()
