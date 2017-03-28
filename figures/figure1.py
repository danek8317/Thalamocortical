import scipy.spatial
import numpy as np
import h5py as h5
from matplotlib import rcParams, cm, ticker, gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
#from KCSD2D import KCSD2D

rcParams.update({'font.size': 8, 'font.family': 'sans-serif'})


def place_electrodes_2D(nx, ny):
    """place nx*ny electrodes next to the column - like an MEA"""
    tot_ele = nx * ny
    zz = np.ones((tot_ele, 1)) * -25.
    xx, yy = np.mgrid[-375:375:np.complex(0, nx), -2050:450:np.complex(0, ny)]
    xx = xx.reshape(tot_ele, 1)
    yy = yy.reshape(tot_ele, 1)
    return np.hstack((xx, yy, zz))


def fetch_mid_pts(h, pop_name):
    """gets the mid points from a file, of a particular population name"""
    all_pts = h['/data/static/morphology/' + pop_name]
    x = (all_pts['x0'] + all_pts['x1']) / 2.
    y = (all_pts['y0'] + all_pts['y1']) / 2.
    z = (all_pts['z0'] + all_pts['z1']) / 2.
    x = x.reshape(x.size, 1)
    y = y.reshape(y.size, 1)
    z = z.reshape(z.size, 1)
    return np.hstack((x, y, z))


def find_map_idx(h, pop_name, field_name):
    """Find the corresponding, locations of morphology in field"""
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
    """computes the inverse distance between src_pos and ele_pos"""
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii, electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(
            src_pos, electrode.reshape(1, 3)).flatten()
    dist_matrix = 1 / dist_matrix  # inverse distance matrix
    return dist_matrix


def pot_vs_time(h, pop_name, field_name, src_pos, ele_pos):
    """returns potentials, at ele_pos, due to src_pos, over time"""
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/' + pop_name + '/' + field_name].value
    src_time = src_time[idx]  # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time) * (1 / (4 * np.pi * 0.3))


def get_all_src_pos(h, pop_names, total_cmpts):
    """Function to compute the positions for a list of populations"""
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        all_srcs[np.sum(total_cmpts[:jj]):np.sum(
            total_cmpts[:jj + 1]), :] = fetch_mid_pts(h, pop_name)
    return all_srcs


def get_extracellular(h, pop_names, time_pts, ele_pos):
    """Fuction to obtain the extracellular potentials at some time points and
    electrode positions"""
    pot_sum = np.zeros((num_ele, time_pts))
    for pop_name in pop_names:
        src_pos = fetch_mid_pts(h, pop_name)
        pot_sum += pot_vs_time(h, pop_name, 'i', src_pos, ele_pos)
        print 'Done extracellular pots for pop_name', pop_name
    return pot_sum


def plot_morp_ele(ax1, src_pos, ele_pos, pot, time_pt):
    """Plots the morphology midpoints and the electrode positions"""
    ax = plt.subplot(121, aspect='equal')
    plt.scatter(src_pos[:, 0], src_pos[:, 1],
                marker='.', alpha=0.7, color='k', s=0.6)
    plt.scatter(ele_pos[:, 0], ele_pos[:, 1],
                marker='x', alpha=0.8, color='r', s=0.9)
    # for tx in range(len(ele_pos[:,0])):
    #    plt.text(ele_pos[tx, 0], ele_pos[tx, 1], str(tx))
    ele_1 = 152
    ele_2 = 148
    plt.scatter(ele_pos[ele_1, 0], ele_pos[ele_1, 1],
                marker='s', color='r', s=14.)
    plt.scatter(ele_pos[ele_2, 0], ele_pos[ele_2, 1],
                marker='s', color='b', s=14.)
    plt.xlabel('X ($\mu$m)')
    plt.ylabel('Y ($\mu$m)')
    plt.title('Morphology and electrodes')
    plt.ylim(ymin=-2150, ymax=550)
    plt.xlim(xmin=-450, xmax=450)

    cbaxes = inset_axes(ax,
                        width="50%",  # width = 10% of parent_bbox width
                        height="17%",  # height : 50%
                        loc=4, borderpad=2.2)

    plt.plot(np.arange(6000), pot[ele_1, :], color='r', linewidth=0.5)
    plt.plot(np.arange(6000), pot[ele_2, :], color='b', linewidth=0.5)

    dummy_line = np.arange(-0.5, 0.5, 0.1)
    plt.plot(np.zeros_like(dummy_line) + time_pt,
             dummy_line, color='black', linewidth=1)

    # ax=plt.gca()
    # ax.arrow(time_pt, -0.1, 0., 0.075, head_width=0.05,
    #          head_length=0.05, width=0.1,
    #          length_includes_head=True, fc='k', ec='k')
    plt.xlim((2750, 3500))  # 4250))
    # plt.xticks(np.arange(3000, 5000, 1000), np.arange(300, 500, 100))
    plt.xticks(np.arange(2750, 3750, 250), np.arange(275, 375, 25))
    plt.ylim((-0.2, 0.12))
    plt.yticks(np.arange(-0.2, 0.1, 0.1), np.arange(-0.2, 0.1, 0.1))
    ax = plt.gca()
    ax.get_yaxis().tick_right()  # set_tick_params(direction='in')

    # inset_plt = plt.plot(cax=cbaxes,
    # cbar = plt.colorbar(cax=cbaxes, ticks=[-lfp_max,0.,lfp_max], orientation='horizontal', format='%.2f')
    # cbar.ax.set_xticklabels([round(-lfp_max,2),str('0 $\mu V$'),round(lfp_max,2)])

    return ax1


def plot_extracellular(ax4, lfp, ele_pos, num_x, num_y, time_pt):
    """Plots the extracellular potentials at a given potentials"""
    ax = plt.subplot(122, aspect='equal')
    lfp *= 1000.
    lfp_max = np.max(np.abs(lfp[:, time_pt]))
    levels = np.linspace(-lfp_max, lfp_max, 16)
    im2 = plt.contourf(ele_pos[:, 0].reshape(num_x, num_y),
                       ele_pos[:, 1].reshape(num_x, num_y),
                       lfp[:, time_pt].reshape(num_x, num_y),
                       levels=levels, cmap=plt.cm.PRGn)

    # cb = plt.colorbar(im2, extend='both')
    # tick_locator = ticker.MaxNLocator(nbins=9, trim=False, prune=None)
    # tick_locator.bin_boundaries(-lfp_max, lfp_max)
    # cb.locator = tick_locator
    # cb.ax.yaxis.set_major_locator(ticker.AutoLocator())
    # cb.update_ticks()
    # cb.ax.set_title('$\mu$V')

    plt.title('Time=' + str(time_pt / 10.) + ' ms')
    plt.xlabel('X ($\mu$m)')
    plt.ylabel('Y ($\mu$m)')
    plt.ylim(ymin=-2150, ymax=550)
    plt.xlim(xmin=-450, xmax=450)
    cbaxes = inset_axes(ax,
                        width="50%",  # width = 10% of parent_bbox width
                        height="2%",  # height : 50%
                        loc=1, borderpad=1.5)
    cbar = plt.colorbar(
        cax=cbaxes, ticks=[-lfp_max, 0., lfp_max], orientation='horizontal',
        format='%.2f')
    cbar.ax.set_xticklabels(
        [round(-lfp_max, 2), str('0 $\mu V$'), round(lfp_max, 2)])
    return ax4


time_pt_interest = 3010
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
cell_range = [0, 1000, 1050, 1140, 1230, 1320,
              1560, 2360, 2560, 3060, 3160, 3260, 3360]
num_cells = np.diff(cell_range) / 10  # 10% MODEL
total_cmpts = list(num_cmpts * num_cells)
pop_names = ['pyrRS23', 'pyrFRB23', 'bask23', 'axax23', 'LTS23',
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
             'bask56', 'axax56', 'LTS56']
h = h5.File('../data/dataset24.h5', 'r')
# h = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/512/12_5/all/traub.h5', 'r')
num_x, num_y = 16, 20
num_ele = num_x * num_y
ele_pos = place_electrodes_2D(num_x, num_y)

time_pts = 6000
pot = get_extracellular(h, pop_names, time_pts, ele_pos)
src_pos = get_all_src_pos(h, pop_names, total_cmpts)
fig = plt.figure()  # figsize=(4,6))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])
ax1 = plt.subplot(gs[:, 0])
ax1 = plot_morp_ele(ax1, src_pos, ele_pos, pot, time_pt_interest)

ax4 = plt.subplot(gs[:, -1])
ax4 = plot_extracellular(ax4, pot, ele_pos, num_x, num_y, time_pt_interest)

# fig, ax3 = plot_csd(fig, ele_pos, pot, time_pt_interest)
plt.tight_layout()
# plt.savefig('fig1.png', dpi=600)
plt.show()
