import numpy as np
import MEAutility as mu
import h5py as h5
from scipy import spatial
import matplotlib.pyplot as plt

# Properties of the Traub's model
cell_range = [0, 1000, 1050, 1140, 1230, 1320,
              1560, 2360, 2560, 3060, 3160, 3260, 3360]
pop_names = ['pyrRS23', 'pyrFRB23', 'bask23', 'axax23', 'LTS23',
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
             'bask56', 'axax56', 'LTS56']
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]


def fetch_simulation(file_loc):
    ''' Properties of the simulated soultions'''
    h = h5.File(file_loc, 'r')
    full_model = True if len(h['/data/static/morphology/pyrRS23']) == 74000 else False
    if full_model:
        num_cells = np.diff(cell_range)  # Full MODEL
    else:
        num_cells = np.diff(cell_range) / 10  # 10% MODEL
    tot_cmpts = list(num_cmpts * num_cells)
    tot_tpts = h['/data/uniform/pyrRS23/i'].shape[1]
    return h, num_cells, tot_cmpts, tot_tpts

def fetch_mid_pts(h, pop_name, test=False):
    ''' Fetch mid_pts from morphology'''
    all_pts = h['/data/static/morphology/' + pop_name]
    x = (all_pts['x0'] + all_pts['x1']) / 2.
    y = (all_pts['y0'] + all_pts['y1']) / 2.
    z = (all_pts['z0'] + all_pts['z1']) / 2.
    x = x.reshape(x.size, 1)
    y = y.reshape(y.size, 1)
    z = z.reshape(z.size, 1)
    locs = np.hstack((z, x, y))
    if test:
        print(pop_name)
        print(locs.shape)
        print('Min x, y, z')
        print(np.min(locs[:, 0]), np.min(locs[:, 1]), np.min(locs[:, 2]))
        print('Max x, y, z')
        print(np.max(locs[:, 0]), np.max(locs[:, 1]), np.max(locs[:, 2]))
    return locs

def fetch_electrodes(test=False, x_offset=0, y_offset=0, z_offset=0):
    ''' fetch locations of electrodes'''
    neuropixels = mu.return_mea('Neuropixels-128')
    ele = neuropixels.positions
    ele[:, 0] += x_offset
    ele[:, 1] += y_offset
    ele[:, 2] += z_offset
    if test:
        print('Electrodes')
        print(ele.shape)
        print('Min x, y, z')
        print(np.min(ele[:, 0]), np.min(ele[:, 1]), np.min(ele[:, 2]))
        print('Max x, y, z')
        print(np.max(ele[:, 0]), np.max(ele[:, 1]), np.max(ele[:, 2]))
    return ele
    
def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/' + pop_name]
    i_data = h['/data/uniform/' + pop_name + '/' + field_name]
    mor_names = mor.dims[0].values()[0]  # entry in /map
    i_names = i_data.dims[0].values()[0]  # entry in /map
    if np.array_equal(i_names, mor_names):
        print('One to one mapping! Neat NSDF!')
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value == entry)[0][0] for entry in mor_names]
    return idx

def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii, electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = spatial.distance.cdist(src_pos,
                                                    electrode.reshape(1, 3)).flatten()
    dist_matrix = 1 / dist_matrix  # inverse distance matrix
    return dist_matrix

def get_all_src_pos(h, total_cmpts):
    """Function to compute the positions for a list of populations"""
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        jjp = jj + 1
        all_srcs[int(np.sum(total_cmpts[:jj])):
                 int(np.sum(total_cmpts[:jjp])), :] = fetch_mid_pts(h, pop_name)
    return all_srcs

def pot_vs_time(h, pop_name, field_name, src_pos, ele_pos):
    '''returns potentials, at ele_pos, due to src_pos, over time, fwd model here'''
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/' + pop_name + '/' + field_name].value
    src_time = src_time[idx]  # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time) * (1 / (4 * np.pi * 0.3))

def make_plot(src_pos, ele_pos, save_filename):
    #fig = plt.figure()
    ax1 = plt.subplot(121, aspect='equal')
    ax1.scatter(src_pos[:, 0], src_pos[:, 1],
                marker='.', alpha=0.7, color='k', s=0.6)
    ax1.scatter(ele_pos[:, 0], ele_pos[:, 1],
                marker='o', color='r', s=0.6)
    ax1.set_xlabel('X ($\mu m)$')
    ax1.set_ylabel('Y ($\mu m)$')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2 = plt.subplot(122, aspect='equal')
    ax2.scatter(src_pos[:, 1], src_pos[:, 2],
                marker='.', alpha=0.7, color='k', s=0.6)
    ax2.scatter(ele_pos[:, 1], ele_pos[:, 2],
                marker='o', color='r', s=0.6)
    ax2.set_xlabel('Y ($\mu m)$')
    ax2.set_ylabel('Z ($\mu m)$')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False) 
    plt.suptitle('Traubs TC model and electrodes')
    plt.savefig(save_filename+'.png', dpi=200)
    plt.show()
    
def save_potentails(file_loc, save_filename, show=False):
    ''' Either show electrodes or compute potentials'''
    ele_pos = fetch_electrodes(x_offset=150)
    num_ele = ele_pos.shape[0]
    h, num_cells, tot_cmpts, tot_tpts = fetch_simulation(file_loc)    
    if show:
        all_src_pos = get_all_src_pos(h, tot_cmpts)
        make_plot(all_src_pos, ele_pos, save_filename)
    else:
        pot_sum = np.zeros((num_ele, tot_tpts))
        for pop_name in pop_names:
            src_pos = fetch_mid_pts(h, pop_name)
            pot_sum += pot_vs_time(h, pop_name, 'i', src_pos, ele_pos)
            print('Done computing extracell pot for pop_name', pop_name)
        # pot_sum = lowpassfilter(pot_sum)
        # an_sign = export_as_neo(ele_pos, pot_sum)
        np.save(save_filename+'.npy', pot_sum)
    h.close()
    

file_loc = '/home/chaitanya/Downloads/pulsestimulus.h5'
save_filename = 'neuropixels_128'

# save_potentails(file_loc, save_filename, show=True)
save_potentails(file_loc, save_filename, show=False)
