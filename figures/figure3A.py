import scipy.spatial
from scipy import signal
import numpy as np
import h5py as h5
import brewer2mpl
from matplotlib import rcParams, cm, gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import sys
sys.path.append('/home/cchintaluri/Thalamocortical/figures/kCSD-python')
from KCSD1D import KCSD1D


rcParams.update({'font.size': 8, 'font.family': 'sans-serif'})
def fetch_mid_pts(h, pop_name):
    all_pts = h['/data/static/morphology/'+pop_name]
    x = (all_pts['x0']+all_pts['x1']) / 2.
    y = (all_pts['y0']+all_pts['y1']) / 2.
    z = (all_pts['z0']+all_pts['z1']) / 2.
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    z = z.reshape(z.size,1)
    return np.hstack((x, y, z))

def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/'+pop_name]
    print pop_name, field_name
    i_data = h['/data/uniform/'+pop_name+'/'+field_name]
    mor_names = mor.dims[0].values()[0] #entry in /map
    i_names = i_data.dims[0].values()[0] #entry in /map
    if np.array_equal(i_names, mor_names):
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value==entry)[0][0] for entry in mor_names]
    return idx

def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii,electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(src_pos, electrode.reshape(1,3)).flatten()
    dist_matrix = 1 / dist_matrix #inverse distance matrix
    return dist_matrix

def pot_vs_time(h, pop_name, field_names, src_pos, ele_pos):
    '''returns potentials, at ele_pos, due to src_pos, over time'''
    ele_src = inv_distance(src_pos, ele_pos).T
    src_time = np.zeros((src_pos.shape[0], 6000))
    for field_name in field_names:
        idx = find_map_idx(h, pop_name, field_name)
        this_field = h['/data/uniform/'+pop_name+'/'+field_name].value
        src_time += this_field[idx] # Order according to correct indices
    return np.dot(ele_src, src_time)*(1 / (4*np.pi*0.3))

def get_all_src_pos(h, pop_names, total_cmpts):
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        all_srcs[np.sum(total_cmpts[:jj]):np.sum(total_cmpts[:jj+1]), :] = fetch_mid_pts(h, pop_name)
    return all_srcs

def get_extracellular(h, pop_names, time_pts, ele_pos, variables):
    pot_sum = np.zeros((num_ele, time_pts))
    for pop_name in pop_names:
        src_pos = fetch_mid_pts(h, pop_name)
        pot_sum += pot_vs_time(h, pop_name, variables, src_pos, ele_pos)
        print 'Done extracellular pots for pop_name, using currents', pop_name, variables
    return pot_sum

def place_electrodes_1D(n):
    '''places n number of electrodes in ID along the column'''
    ele_x = np.ones((n,1))*25.
    ele_y = np.linspace(-2050., 450.0, num=n).reshape(n,1)
    ele_z = np.ones((n,1))*25.
    return np.hstack((ele_x, ele_y, ele_z))

def lowpassfilter(data):
    sampling_rate = 10000 #10kHz, corresp to 0.1ms
    nyq = 0.5 * sampling_rate
    normal_cutoff = 100 / nyq
    b, a = signal.butter(2, normal_cutoff, 'low', analog=False)
    return signal.lfilter(b, a, data)

def plot_potential(ax, lfp, max_val, title):
    norm = cm.colors.Normalize(vmax=max_val, vmin=-max_val, clip=False)
    im = plt.imshow(lfp[::-1], aspect='auto', norm=norm, interpolation='nearest', cmap=plt.cm.PRGn)
    plt.xlim((2750, 4250))
    plt.xticks(np.arange(3000, 5000, 1000), np.arange(300, 500, 100))
    #plt.ylabel('Electrode depth ($\mu$m)')
    #plt.xlabel('Time (ms)')
    plt.title(title)#, fontweight="bold", fontsize=12)
    #plt.yticks(np.arange(28))
    plt.gca().set_yticks(np.arange(28))
    plt.gca().set_yticklabels(np.arange(1,29))
    # for label in plt.gca().yaxis.get_ticklabels()[::2]:
    #     label.set_visible(False)
    #plt.xlim(xmin=2500, xmax=4500)
    #plt.colorbar(extend='both')
    cbaxes = inset_axes(ax,
                        width="40%",  # width = 10% of parent_bbox width
                        height="3%",  # height : 50%
                        loc=1, borderpad=1 )
    cbar = plt.colorbar(cax=cbaxes, ticks=[-max_val,0.,max_val], orientation='horizontal', format='%.2f')
    cbar.ax.set_xticklabels([round(-max_val,2),str('0 $\mu V$'),round(max_val,2)])
    return ax, im



def plot_raster(h, ax, title):
    set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
    set2.insert(0, 'crimson')
    set2.append('dodgerblue')
    cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360,3460,3560]
    pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
                 'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
                 'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']
    inh_cells = ['bask23','axax23','LTS23', 
                 'bask56', 'axax56', 'LTS56', 'nRT']
    cells = np.array(cell_range)/ 10 #10% MODEL
    cells_spike = np.diff(cells)
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
            #ax.plot(x, y, '.', color=color)
            ax.scatter(x, y, s=7, alpha=0.8,facecolor=color, 
                       marker=marker, linewidth='0.2', edgecolor='gray') #edgecolor='none'
            count -= intra_cell_gap
        count -= intra_pop_gap
    mid_pts = cells_spike / 2
    tick_pts = -1*(mid_pts + cells[:-1])*intra_cell_gap
    for ii in range(len(pop_names)):
        tick_pts[ii] -= (intra_pop_gap*ii) 
    ax.set_yticks(tick_pts)
    ax.set_yticklabels(pop_names)
    plt.plot(np.zeros((5000))+300, np.arange(-4500,500), color='k')
    plt.xlim((275, 425))
    plt.ylim((-2400 ,50))
    plt.xticks(np.arange(300, 500, 100), np.arange(300, 500, 100))
    plt.title(title)#, fontweight="bold", fontsize=12)
    return ax

def get_specific(ax, h, pop_names, time_pts, ele_pos, variables, title):
    if title == 'Other':
        pot_all_but = get_extracellular(h, pop_names, time_pts, ele_pos, variables)
        pot_all = get_extracellular(h, pop_names, time_pts, ele_pos, ['i'])
        pot = pot_all - pot_all_but
    else:
        pot = get_extracellular(h, pop_names, time_pts, ele_pos, variables)
    lfp = lowpassfilter(pot)
    lfp_max = np.max(np.abs(lfp[:, 2750:4250]))
    ax,im = plot_potential(ax, lfp, lfp_max, title)
    return ax, im, lfp

def set_axis(ax, letter=None):
    if letter is not None:
        ax.text(0.05, 1.025, letter, fontsize=12, weight='bold', transform=ax.transAxes)
    #ax.set_aspect('equal')
    return ax

def add_second_yaxis(ax):
    #ax.autoscale(False)
    ax2 = ax.twinx()
    ax2.set_ylabel("Electrode position (mm)")
    ele_pos = place_electrodes_1D(28)
    ax2.set_yticks(np.arange(28))
    ax2.set_yticklabels(np.round(ele_pos[:,1]/1000.,2))
    for label in ax2.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    return


def place_electrodes_2D(nx, ny):
    """place nx*ny electrodes next to the column - like an MEA"""
    tot_ele = nx*ny
    zz = np.ones((tot_ele,1))*-25.
    xx, yy = np.mgrid[-375:375:np.complex(0,nx), -2050:450:np.complex(0,ny)]
    xx = xx.reshape(tot_ele, 1)
    yy = yy.reshape(tot_ele, 1)
    return np.hstack((xx, yy, zz))

def fetch_mid_pts(h, pop_name):
    """gets the mid points from a file, of a particular population name"""
    all_pts = h['/data/static/morphology/'+pop_name]
    x = (all_pts['x0']+all_pts['x1']) / 2.
    y = (all_pts['y0']+all_pts['y1']) / 2.
    z = (all_pts['z0']+all_pts['z1']) / 2.
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    z = z.reshape(z.size,1)
    return np.hstack((x, y, z))

def find_map_idx(h, pop_name, field_name):
    """Find the corresponding, locations of morphology in field"""
    mor = h['/data/static/morphology/'+pop_name]
    i_data = h['/data/uniform/'+pop_name+'/'+field_name]
    mor_names = mor.dims[0].values()[0] #entry in /map
    i_names = i_data.dims[0].values()[0] #entry in /map
    if np.array_equal(i_names, mor_names):
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value==entry)[0][0] for entry in mor_names]
    return idx

def inv_distance(src_pos, ele_pos):
    """computes the inverse distance between src_pos and ele_pos"""
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii,electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(src_pos, electrode.reshape(1,3)).flatten()
    dist_matrix = 1 / dist_matrix #inverse distance matrix
    return dist_matrix

def pot_vs_time_2D(h, pop_name, field_name, src_pos, ele_pos):
    """returns potentials, at ele_pos, due to src_pos, over time"""
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/'+pop_name+'/'+field_name].value
    src_time = src_time[idx] # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time)*(1 / (4*np.pi*0.3))

def get_all_src_pos(h, pop_names, total_cmpts):
    """Function to compute the positions for a list of populations"""
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        all_srcs[np.sum(total_cmpts[:jj]):np.sum(total_cmpts[:jj+1]), :] = fetch_mid_pts(h, pop_name)
    return all_srcs

def get_extracellular_2D(h, pop_names, time_pts, ele_pos):
    """Fuction to obtain the extracellular potentails at some time points and electrode positions"""
    num_ele = ele_pos.shape[0]
    pot_sum = np.zeros((num_ele, time_pts))
    for pop_name in pop_names:
        src_pos = fetch_mid_pts(h, pop_name)
        pot_sum += pot_vs_time_2D(h, pop_name, 'i', src_pos, ele_pos)
        print 'Done extracellular pots for pop_name', pop_name
    return pot_sum

def plot_morp_ele(ax, src_pos, ele_pos, pot, time_pt):
    """Plots the morphology midpoints and the electrode positions""" 
    #ax = plt.subplot(121, aspect='equal')
    plt.scatter(src_pos[:, 0], src_pos[:, 1], marker='.', alpha=0.7, color='k', s=0.6)
    plt.scatter(ele_pos[:, 0], ele_pos[:, 1], marker='x', alpha=0.8, color='r', s=0.9)
    #for tx in range(len(ele_pos[:,0])):
    #    plt.text(ele_pos[tx, 0], ele_pos[tx, 1], str(tx))
    ele_1 = 152
    ele_2 = 148
    plt.scatter(ele_pos[ele_1, 0], ele_pos[ele_1, 1], marker='s', color='r', s=14.)
    plt.scatter(ele_pos[ele_2, 0], ele_pos[ele_2, 1], marker='s', color='b', s=14.)
    plt.xlabel('X ($\mu$m)')
    plt.ylabel('Y ($\mu$m)')
    plt.title('Morphology, electrodes')
    plt.ylim(ymin=-2150,ymax=550)
    plt.xlim(xmin=-450,xmax=450)

    cbaxes = inset_axes(ax,
                        width="50%",  # width = 10% of parent_bbox width
                        height="17%",  # height : 50%
                        loc=4, borderpad=2.2)
    

    plt.plot(np.arange(6000), pot[ele_1, :], color='r', linewidth=0.5)
    plt.plot(np.arange(6000), pot[ele_2, :], color='b', linewidth=0.5)

    dummy_line = np.arange(-0.5, 0.5, 0.1)
    plt.plot(np.zeros_like(dummy_line)+time_pt, dummy_line, color='black', linewidth=1) 

    # ax=plt.gca()
    # ax.arrow(time_pt, -0.1, 0., 0.075, head_width=0.05,
    #          head_length=0.05, width=0.1,
    #          length_includes_head=True, fc='k', ec='k')
    plt.xlim((2750, 3500)) #4250))
    #plt.xticks(np.arange(3000, 5000, 1000), np.arange(300, 500, 100))
    plt.xticks(np.arange(2750, 3750, 250), np.arange(275, 375, 25))
    plt.ylim((-0.2, 0.12))
    plt.yticks(np.arange(-0.2, 0.1, 0.1),np.arange(-0.2, 0.1, 0.1))
    ax = plt.gca()
    ax.get_yaxis().tick_right()#set_tick_params(direction='in')


    # inset_plt = plt.plot(cax=cbaxes, 
    # cbar = plt.colorbar(cax=cbaxes, ticks=[-lfp_max,0.,lfp_max], orientation='horizontal', format='%.2f')
    # cbar.ax.set_xticklabels([round(-lfp_max,2),str('0 $\mu V$'),round(lfp_max,2)])

    return ax

def plot_extracellular(ax, lfp, ele_pos, num_x, num_y, time_pt):
    """Plots the extracellular potentials at a given potentials"""

    lfp *= 1000.
    lfp_max = np.max(np.abs(lfp[:, time_pt]))
    levels = np.linspace(-lfp_max, lfp_max, 16)
    im2 = plt.contourf(ele_pos[:,0].reshape(num_x, num_y), 
                       ele_pos[:,1].reshape(num_x, num_y), 
                       lfp[:,time_pt].reshape(num_x,num_y), 
                       levels=levels, cmap=plt.cm.PRGn)


    # cb = plt.colorbar(im2, extend='both')
    # tick_locator = ticker.MaxNLocator(nbins=9, trim=False, prune=None)
    # #tick_locator.bin_boundaries(-lfp_max, lfp_max)
    # cb.locator = tick_locator
    # #cb.ax.yaxis.set_major_locator(ticker.AutoLocator())
    # cb.update_ticks()
    # cb.ax.set_title('$\mu$V')

    plt.title('Time='+str(time_pt/10.)+' ms')
    plt.xlabel('X ($\mu$m)')
    plt.ylabel('Y ($\mu$m)')
    plt.ylim(ymin=-2150,ymax=550)
    plt.xlim(xmin=-450,xmax=450)
    cbaxes = inset_axes(ax,
                        width="50%",  # width = 10% of parent_bbox width
                        height="2%",  # height : 50%
                        loc=1, borderpad=1.5)
    cbar = plt.colorbar(cax=cbaxes, ticks=[-lfp_max,0.,lfp_max], orientation='horizontal', format='%.2f')
    cbar.ax.set_xticklabels([round(-lfp_max,2),str('0 $\mu V$'),round(lfp_max,2)])
    return ax


def plot_csd(ax, lfp, max_val, title, start_t, end_t):
    norm = cm.colors.Normalize(vmax=max_val, vmin=-max_val, clip=False)
    im = plt.imshow(lfp[::-1], aspect='auto', norm=norm, interpolation='nearest', cmap=plt.cm.bwr_r)
    #plt.xlim((2750, 4250))
    plt.xticks(np.arange(3000-start_t, 5000-start_t, 1000), np.arange(300, 500, 100))
    #plt.ylabel('Electrode depth ($\mu$m)')
    plt.xlabel('Time (ms)')
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.title(title)#, fontweight="bold", fontsize=12)
    #plt.yticks(np.arange(28))
    # plt.gca().set_yticks(np.arange(28))
    # plt.gca().set_yticklabels(np.arange(1,29))
    for label in plt.gca().yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    #plt.xlim(xmin=2500, xmax=4500)
    #plt.colorbar(extend='both')
    cbaxes = inset_axes(ax,
                        width="40%",  # width = 10% of parent_bbox width
                        height="3%",  # height : 50%
                        loc=1, borderpad=1 )
    cbar = plt.colorbar(cax=cbaxes, ticks=[-max_val,0.,max_val], orientation='horizontal', format='%.2f')
    cbar.ax.set_xticklabels(['source',str('0'),'sink'])
    return ax, im


def get_csd(ax, ele_pos, lfp):
    ele_pos = ele_pos[:, 1][:, np.newaxis]/1000.
    start_t = 2750
    end_t = 4250
    # k = KCSD1D(ele_pos, lfp[:, 3200:3210]*1000.)
    # k.cross_validate(Rs=np.array([0.05, 0.1, 0.41, 0.8, 1.6]), lambdas=np.logspace(15,-25, 35))
    k = KCSD1D(ele_pos, lfp[:, start_t:end_t]*1000., 
               R_init=0.41, lambd=1e-05)
    #k.cross_validate(Rs=np.array([0.41]), lambdas=np.logspace(15,-25, 35))
    est_csd = k.values()
    csd_max = np.max(np.abs(est_csd))
    ax, im = plot_csd(ax, est_csd, csd_max, 'CSD', start_t, end_t)
    return ax, im

def consolidated_figure(h, pop_names, ele_pos, time_pts, fig, variable_dict):
    z_steps = 2 #4 for the plots 
    height_ratios = [1 for i in range(z_steps)]
    width_ratios = [0.5, 0.5, 1.2]
    #height_ratios.append(0.07) #height of colorbar
    gs = gridspec.GridSpec(z_steps, 3, height_ratios=height_ratios, width_ratios=width_ratios)
    #Morphology
    ax = plt.subplot(gs[0, 0])
    time_pt_interest = 3010
    num_x_2d, num_y_2d = 16, 20
    num_ele = num_x_2d * num_y_2d
    ele_pos_2D = place_electrodes_2D(num_x_2d, num_y_2d)
    time_pts = 6000
    pot_2D = get_extracellular_2D(h, pop_names, time_pts, ele_pos_2D)
    src_pos = get_all_src_pos(h, pop_names, total_cmpts)
    morph_plot = plot_morp_ele(ax, src_pos, ele_pos_2D, pot_2D, time_pt_interest)
    set_axis(ax, letter='A')
    #2D potentails
    ax = plt.subplot(gs[0, 1])    
    pot_2d_ax = plot_extracellular(ax, pot_2D, ele_pos_2D, num_x_2d, num_y_2d, time_pt_interest)
    set_axis(ax, letter='B')
    #Spikes
    ax = plt.subplot(gs[0, 2])
    spike_ax = plot_raster(h, ax, 'Spikes')
    set_axis(ax, letter='C')
    ##LFP
    ax = plt.subplot(gs[1, 0:-1])
    key_val = 'LFP'
    ax, im, lfp = get_specific(ax, h, pop_names, time_pts, ele_pos, variable_dict[key_val], key_val)
    set_axis(ax, letter='D')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Electrode number')
    # CSD
    ax = plt.subplot(gs[1, 2])
    ax, im = get_csd(ax, ele_pos, lfp)
    set_axis(ax, letter='E')
    ax.set_xlabel('Time (ms)')
    

    return fig
    
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360]
num_cells = np.diff(cell_range) / 10 #10% MODEL
total_cmpts = list(num_cmpts*num_cells)
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56']
h = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/1/traub_syn.h5', 'r')

# #IGNORED MAX_VAL_DICT
# max_val_dict = {'i':0.05, 'i_gaba_a': 0.05, 'i_cap':0.05, 
#                 'i_pas':0.05, 'i_k':0.5, 'i_na':0.5, 'i_ca':0.5, 
#                 'i_cat':0.05, 'i_ar':0.05, 'AMPA+NMDA':0.5, 'depol_curr':0.05}

# title_dict = {'i':'LFP', 'i_gaba_a': 'GABA A', 'i_cap':'Capacitive', 
#               'i_pas':'Passive', 'i_k':'Potassium', 'i_na':'Sodium', 
#               'i_ca':'Calcium', 'i_cat':'Calcium T type', 'i_ar':'Anomalous rectifier', 
#               'AMPA+NMDA':'AMPA+NMDA', 'depol_curr':'Depolarizing Current'}

variable_dict = {'LFP':['i'], 'AMPA+NMDA':['i_AMPA', 'i_NMDA'], 'GABA A':['i_gaba_a'],
                 'Capacitive':['i_cap'], 'Potassium':['i_k'], 'Passive':['i_pas'],
                 'Sodium':['i_na'], 'Calcium T type':['i_cat', 'i_cat_a'], 
                 'Anomalous rectifier':['i_ar'], 'Calcium':['i_ca'], 
                 'Other':['i_AMPA', 'i_NMDA', 'i_gaba_a', 'i_cap', 'i_k', 'i_pas', 'i_na', 'i_cat', 'i_cat_a', 'i_ar', 'i_ca']
                 }

num_ele = 28
ele_pos = place_electrodes_1D(num_ele)
time_pts = 6000

#pot = get_extracellular(h, pop_names, time_pts, ele_pos, 'i')
#lfp = lowpassfilter(pot)
#fig = plt.figure()#figsize=(4,6))
#ax = subplot(111)
#ax1, im = plot_potential(ax, lfp, 0.05, 'LFP')
#ax1 = plot_raster(h, ax, 'Spikes')

fig = plt.figure(figsize=(12,12))
fig = consolidated_figure(h, pop_names, ele_pos, time_pts, fig, variable_dict)
plt.tight_layout()
plt.savefig('fig3A.png', dpi=300)
#plt.show()
