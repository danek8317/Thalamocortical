import scipy.spatial
from scipy import signal
import numpy as np
import h5py as h5
import brewer2mpl
from matplotlib import rcParams, cm  # , gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

rcParams.update({'font.size': 8, 'font.family': 'sans-serif'})


def fetch_mid_pts(h, pop_name):
    all_pts = h['/data/static/morphology/' + pop_name]
    x = (all_pts['x0'] + all_pts['x1']) / 2.
    y = (all_pts['y0'] + all_pts['y1']) / 2.
    z = (all_pts['z0'] + all_pts['z1']) / 2.
    x = x.reshape(x.size, 1)
    y = y.reshape(y.size, 1)
    z = z.reshape(z.size, 1)
    return np.hstack((x, y, z))


def fetch_soma_idx(h, pop_name, cell_numbers):
    ''' cell_numbers is a list of int of cells'''
    mor = h['/data/static/morphology/' + pop_name]
    mor_names = mor.dims[0].values()[0]  # entry in /map
    idxs = [0] * len(cell_numbers)
    for jj, val in enumerate(mor_names):
        cell, cmpt = val.split('/')
        cell = int(cell)
        if (cell in cell_numbers) and cmpt == '1':
            idxs[cell_numbers.index(cell)] = jj
    return idxs


def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/' + pop_name]
    i_data = h['/data/uniform/' + pop_name + '/' + field_name]
    mor_names = mor.dims[0].values()[0]  # entry in /map
    i_names = i_data.dims[0].values()[0]  # entry in /map
    if np.array_equal(i_names, mor_names):
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value == entry)[0][0] for entry in mor_names]
    return idx


def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii, electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(
            src_pos, electrode.reshape(1, 3)).flatten()
    dist_matrix = 1 / dist_matrix  # inverse distance matrix
    return dist_matrix


def pot_vs_time(h, pop_name, field_name, src_pos, ele_pos):
    '''returns potentials, at ele_pos, due to src_pos, over time'''
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/' + pop_name + '/' + field_name].value
    src_time = src_time[idx]  # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time) * (1 / (4 * np.pi * 0.3))


def get_all_src_pos(h, pop_names, total_cmpts):
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        all_srcs[np.sum(total_cmpts[:jj]):np.sum(
            total_cmpts[:jj + 1]), :] = fetch_mid_pts(h, pop_name)
    return all_srcs


def get_extracellular(h, pop_names, time_pts, ele_pos, variable):
    pot_sum = np.zeros((num_ele, time_pts))
    for pop_name in pop_names:
        src_pos = fetch_mid_pts(h, pop_name)
        pot_sum += pot_vs_time(h, pop_name, variable, src_pos, ele_pos)
        print('Done extracellular pots for pop_name, using currents', pop_name, variable)
    return pot_sum


def place_electrodes_1D(n):
    '''places n number of electrodes in ID along the column'''
    ele_x = np.ones((n, 1)) * 25.
    ele_y = np.linspace(-2050., 450.0, num=n).reshape(n, 1)
    ele_z = np.ones((n, 1)) * 25.
    return np.hstack((ele_x, ele_y, ele_z))


def lowpassfilter(data):
    sampling_rate = 10000  # 10kHz, corresp to 0.1ms
    nyq = 0.5 * sampling_rate
    normal_cutoff = 100 / nyq
    b, a = signal.butter(2, normal_cutoff, 'low', analog=False)
    return signal.lfilter(b, a, data)


def plot_potential(ax, lfp, max_val, title):
    norm = cm.colors.Normalize(vmax=max_val, vmin=-max_val, clip=False)
    im = plt.imshow(
        lfp[::-1],
        aspect='auto',
        norm=norm,
        interpolation='nearest',
        cmap=plt.cm.PRGn)

    # plt.xticks(np.arange(3000, 5000, 1000), np.arange(300, 500, 100))
    # plt.ylabel('Electrode depth ($\mu$m)')
    # plt.xlabel('Time (ms)')
    plt.title(title, fontweight="bold", fontsize=12)
    # plt.xlim(xmin=2500, xmax=4500)
    # plt.colorbar(extend='both')
    plt.gca().set_yticks(np.arange(28))
    plt.gca().set_yticklabels(np.arange(1, 29))
    for label in plt.gca().yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    cbaxes = inset_axes(ax,
                        width="40%",  # width = 10% of parent_bbox width
                        height="3%",  # height : 50%
                        loc=1, borderpad=1)
    cbar = plt.colorbar(
        cax=cbaxes,
        ticks=[-max_val,
               0.,
               max_val],
        orientation='horizontal',
        format='%.2f')
    cbar.ax.set_xticklabels(
        [round(-max_val, 2), str('0 $\mu V$'), round(max_val, 2)])
    return ax, im


def get_specific(ax, h, pop_names, time_pts, ele_pos, variable, title):
    if title == 'LFP':
        pot = get_extracellular(
            h,
            pop_names,
            time_pts,
            ele_pos,
            variable)  # LOWPASS
        lfp = lowpassfilter(pot)
    else:
        lfp = get_extracellular(h, pop_names, time_pts, ele_pos, variable)
    ax, im = plot_potential(ax, lfp, 0.05, title)
    return ax, im


def set_axis(ax, letter=None):
    if letter is not None:
        ax.text(
            0.05,
            1.025,
            letter,
            fontsize=12,
            weight='bold',
            transform=ax.transAxes)
    # ax.set_aspect('equal')
    return ax


def add_second_yaxis(ax):
    # ax.autoscale(False)
    ax2 = ax.twinx()
    ax2.set_ylabel("Electrode position (um)")
    ele_pos = place_electrodes_1D(28)
    ax2.set_yticks(np.arange(28))
    ax2.set_yticklabels(np.round(ele_pos[:, 1], 2))
    for label in ax2.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    return


def figure4A(h_dict, pop_names, ele_pos, time_pts, fig):
    ax = plt.subplot(121)
    title = 'Extracellular'
    h = h_dict[title]

    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i', title)
    set_axis(ax, 'A')
    ax.set_ylabel('Electrode number')
    ax.set_xlabel('Time (0.1 ms)')

    add_second_yaxis(ax)
    cell_range = [
        0,
        1000,
        1050,
        1140,
        1230,
        1320,
        1560,
        2360,
        2560,
        3060,
        3160,
        3260,
        3360,
        3460,
        3560]

    # assign color dict
    set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
    set2.insert(0, 'crimson')
    set2.append('dodgerblue')
    color_pop_names = ['pyrRS23', 'pyrFRB23', 'bask23', 'axax23', 'LTS23',
                       'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
                       'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']
    color_dict = {}
    for cc, color_pop_name in enumerate(color_pop_names):
        color_dict[color_pop_name] = set2[cc]

    inh_cells = ['bask23', 'axax23', 'LTS23',
                 'bask56', 'axax56', 'LTS56', 'nRT']
    cells = np.array(cell_range) / 10  # 10% MODEL
    cells_spike = np.diff(cells)

    ax2 = plt.subplot(122)

    for ii, pop_name in enumerate(pop_names):
        cell_count = cells_spike[ii]
        spike_array = h['/data/event/' + pop_name + '/spikes']
        all_spikes = spike_array.value
        cell_names = list(spike_array.dims[0].values()[0].value)
        all_cmpts = fetch_mid_pts(h, pop_name)
        soma_idx = fetch_soma_idx(h, pop_name, cell_names)

        soma_x = all_cmpts[soma_idx, 0]
        soma_y = all_cmpts[soma_idx, 1]
        soma_z = all_cmpts[soma_idx, 2]

        color = color_dict[pop_name]
        if pop_name in inh_cells:
            marker = 'v'
        else:
            marker = '^'

        # distances between somas and electrode
        dists = np.sqrt((soma_x - 25.) ** 2 + (soma_z - 25.) ** 2)
        dists_diff = max(dists) - min(dists)
        dists = (dists - min(dists)) / dists_diff

        for kk in range(all_spikes.shape[0]):
            xx = all_spikes[kk]
            yy = np.zeros_like(xx) + soma_y[kk]
            ax2.scatter(xx * 10, yy, s=15, alpha=dists[kk], facecolor=color,
                        marker=marker, linewidth='0.2', edgecolor='black')  # edgecolor='none'

    # ele_pos = place_electrodes_1D(28)
    # ax2.set_yticks(ele_pos[:,1])
    # ax2.set_yticklabels(np.round(ele_pos[:,1]/1000.,2))
    plt.title('Spikes', fontweight="bold", fontsize=12)
    ax2.set_ylim((-2050., 450.0))
    # ax2.set_yticks([])
    plt.xlabel('Time (ms)')
    # for label in ax2.yaxis.get_ticklabels()[::2]:
    #     label.set_visible(False)

    ax.set_xlim((3000, 3500))
    ax2.set_xlim((3000, 3500))
    ax.xaxis.grid(True)
    ax2.xaxis.grid(True)

    return fig

# def plot_raster(h, ax, title):
#     set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
#     set2.insert(0, 'crimson')
#     set2.append('dodgerblue')
#     cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360,3460,3560]
#     pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23',
#                  'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
#                  'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']
#     inh_cells = ['bask23','axax23','LTS23',
#                  'bask56', 'axax56', 'LTS56', 'nRT']
# cells = np.array(cell_range)/ 10 #10% MODEL
#     cells_spike = np.diff(cells)
#     intra_cell_gap = 6
#     intra_pop_gap = 15
#     count = 0
#     for ii,pop_name in enumerate(pop_names):
#         cell_count = cells_spike[ii]
#         all_spikes = h['/data/event/'+pop_name+'/spikes'].value
#         int_spikes = all_spikes[0:cell_count]
#         color = set2[ii]
#         if pop_name in inh_cells:
#             marker = 'v'
#         else:
#             marker = '^'
#         for kk in range(cell_count):
#             x = int_spikes[kk]
#             y = np.zeros_like(x) + count
# ax.plot(x, y, '.', color=color)
#             ax.scatter(x, y, s=7, alpha=0.8,facecolor=color,
# marker=marker, linewidth='0.2', edgecolor='gray') #edgecolor='none'
#             count -= intra_cell_gap
#         count -= intra_pop_gap
#     mid_pts = cells_spike / 2
#     tick_pts = -1*(mid_pts + cells[:-1])*intra_cell_gap
#     for ii in range(len(pop_names)):
#         tick_pts[ii] -= (intra_pop_gap*ii)
#     ax.set_yticks(tick_pts)
#     ax.set_yticklabels(pop_names)
#     plt.plot(np.zeros((5000))+300, np.arange(-4500,500), color='k')
#     plt.xlim((275, 425))
#     plt.ylim((-2400 ,50))
#     plt.xticks(np.arange(300, 500, 100), np.arange(300, 500, 100))
# plt.title(title, fontweight="bold", fontsize=12)
#     return ax


# def consolidated_figure(h_dict, pop_names, ele_pos, time_pts, fig):
# z_steps = 2 #4 for the plots 1 for colorbar
#     height_ratios = [1 for i in range(z_steps)]
# height_ratios.append(0.07) #height of colorbar
#     gs = gridspec.GridSpec(z_steps, 3, height_ratios=height_ratios)
# extracellular
#     ax = plt.subplot(gs[0, 0])
#     title = 'Extracellular'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     set_axis(ax, 'A')
#     ax.set_ylabel('Electrode number')
# LFP
#     ax = plt.subplot(gs[0, 1])
#     title = 'LFP'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     set_axis(ax, 'B')
# Only Passive
#     ax = plt.subplot(gs[0, 2])
#     title = 'Passive only'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     add_second_yaxis(ax)
#     set_axis(ax, 'C')
# Only Passive soma and axon
#     ax = plt.subplot(gs[1, 0])
#     title = 'Passive soma and axon'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     set_axis(ax, 'D')
#     ax.set_xlabel('Time (ms)')
#     ax.set_ylabel('Electrode number')
# cax = plt.subplot(gs[2, 0])
# cbar = plt.colorbar(im, cax=cax, orientation='horizontal', extend='both')
# passive axon
#     ax = plt.subplot(gs[1, 1])
#     title = 'Passive axon'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     set_axis(ax, 'E')
#     ax.set_xlabel('Time (ms)')
# cax2 = plt.subplot(gs[2, 1])
# cbar2 = plt.colorbar(im, cax=cax2, orientation='horizontal', extend='both')
# closed fast Na
#     ax = plt.subplot(gs[1, 2])
#     title = 'Closed fast Sodium'
#     h = h_dict[title]
#     ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
#     set_axis(ax, 'F')
#     ax.set_xlabel('Time (ms)')
#     add_second_yaxis(ax)
# cax3 = plt.subplot(gs[2, 2])
# cbar3 = plt.colorbar(im, cax=cax3, orientation='horizontal', extend='both')
#     return fig

num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
cell_range = [
    0,
    1000,
    1050,
    1140,
    1230,
    1320,
    1560,
    2360,
    2560,
    3060,
    3160,
    3260,
    3360]
num_cells = np.diff(cell_range) / 10  # 10% MODEL
total_cmpts = list(num_cmpts * num_cells)

# pop_names = ['pyrRS23','pyrFRB23']#,'bask23','axax23','LTS23']

# pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23',
#              'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
#              'bask56', 'axax56', 'LTS56']

# pop_names = ['tuftIB5', 'tuftRS5', 'nontuftRS6']
pop_names = ['nontuftRS6']  # , 'nontuftRS6']
# pop_names = ['nontuftRS6']
# pop_names = ['pyrFRB23']

h24 = h5.File('../data/dataset24.h5', 'r')
h25 = h5.File('../data/dataset25.h5', 'r')
h26 = h5.File('../data/dataset26.h5', 'r')
h27 = h5.File('../data/dataset27.h5', 'r')
h28 = h5.File('../data/dataset28.h5', 'r')
hs = [h24, h25, h26, h27, h28]
h_dict = {'Extracellular': h24,
          'LFP': h24,
          'Passive only': h25,
          'Passive soma and axon': h26,
          'Passive axon': h27,
          'Closed fast Sodium': h28}


max_val_dict = {'i': 0.05, 'i_gaba_a': 0.05, 'i_cap': 0.05,
                'i_pas': 0.05, 'i_k': 0.5, 'i_na': 0.5, 'i_ca': 0.5,
                'i_cat': 0.05, 'i_ar': 0.05, 'AMPA+NMDA': 0.5}
title_dict = {'i': 'LFP', 'i_gaba_a': 'GABA A', 'i_cap': 'Capacitive',
              'i_pas': 'Passive', 'i_k': 'Potassium', 'i_na': 'Sodium',
              'i_ca': 'Calcium', 'i_cat': 'CAT', 'i_ar': 'AR', 'AMPA+NMDA': 'AMPA+NMDA'}
num_ele = 28
ele_pos = place_electrodes_1D(num_ele)
time_pts = 6000
fig = plt.figure(figsize=(16, 8))
# fig = consolidated_figure(h_dict, pop_names, ele_pos, time_pts, fig)
fig = figure4A(h_dict, pop_names, ele_pos, time_pts, fig)
[h.close() for h in hs]  # close all files memory overflow
plt.tight_layout()
plt.savefig('fig4A_nontuftRS6.png', dpi=600)
# plt.savefig('fig4A_nontuftRS6.png', dpi=600)
# plt.show()
