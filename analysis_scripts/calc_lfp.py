from scipy import signal
import numpy as np
import h5py as h5
import scipy.spatial
from lfp_parameters import *


def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/' + pop_name]
    i_data = h['/data/uniform/' + pop_name + '/' + field_name]
    mor_names = mor.dims[0].values()[0]  # entry in /map
    i_names = i_data.dims[0].values()[0]  # entry in /map
    if np.array_equal(i_names, mor_names):
        print 'One to one mapping! Neat NSDF!'
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value == entry)[0][0] for entry in mor_names]
    return idx


def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii, electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(src_pos,
                                                          electrode.reshape(1, 3)).flatten()
    dist_matrix = 1 / dist_matrix  # inverse distance matrix
    return dist_matrix


def pot_vs_time(h, pop_name, field_name, src_pos, ele_pos):
    '''returns potentials, at ele_pos, due to src_pos, over time'''
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/' + pop_name + '/' + field_name].value
    src_time = src_time[idx]  # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time) * (1 / (4 * np.pi * 0.3))


def lowpassfilter(data):
    '''Lowpass filter for the extracellular potentails 100Hz'''
    sampling_rate = 10000  # 10kHz, corresp to 0.1ms
    nyq = 0.5 * sampling_rate
    normal_cutoff = 100 / nyq
    b, a = signal.butter(2, normal_cutoff, 'low', analog=False)
    return signal.lfilter(b, a, data)


def export_as_neo(ele_pos, pots):
    '''This fuction exports the potentials recorded into neo object list
    this can then be used in third party software such as elephant

    elephant provides rich additions such as spike sorting, csd
    analysis etc.'''
    import neo
    import quantities as pq
    an_sigs = []
    for ii in range(len(pots)):
        rc = neo.RecordingChannel()
        rc.coordinate = ele_pos[ii] * pq.um
        asig = neo.AnalogSignal(pots[ii] * pq.mV, sampling_rate=10000 * pq.Hz)
        rc.analogsignals = [asig]
        rc.create_relationship()
        an_sigs.append(asig)
    print 'Finished converting to neo object list'
    return an_sigs

pot_sum = np.zeros((num_ele, 6000))
for pop_name in pop_names:
    src_pos = fetch_mid_pts(h, pop_name)
    pot_sum += pot_vs_time(h, pop_name, 'i', src_pos, ele_pos)
    print 'Done for pop_name', pop_name
# pot_sum = lowpassfilter(pot_sum)
# an_sign = export_as_neo(ele_pos, pot_sum)
np.save(pot_filename, pot_sum)
h.close()
